/*
    LoudnessMeasurer
    Copyright (c) 2014 Ricci Adams

    Heavily based on libebur128
    Copyright (c) 2011 Jan Kokemüller

    Permission is hereby granted, free of charge, to any person obtaining a copy of
    this software and associated documentation files (the "Software"), to deal in
    the Software without restriction, including without limitation the rights to
    use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
    the Software, and to permit persons to whom the Software is furnished to do so,
    subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
    FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
    COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
    IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "LoudnessMeasurer.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#if __x86__
#include <xmmintrin.h>
#endif

#include <Accelerate/Accelerate.h>
#include <AudioToolbox/AudioToolbox.h>

typedef struct LoudnessMeasurerChannel {
    float  *_bufferPre;
    double *_bufferPost;
    size_t  _bufferFrames;
    size_t  _bufferIndex;

    double *_scratch;
    double *_scratch2;

    double _filterState[2][4];

    double _samplePeak;

    double *_blocks;
    size_t  _blocksCapacity;
    size_t  _blocksCount;

    float  *_overview;
    size_t  _overviewCapacity;
    size_t  _overviewCount;

    // How many frames are needed for a gating block. Will correspond to 400ms
    // of audio at initialization, and 100ms after the first block (75% overlap
    // as specified in the 2011 revision of BS1770). */
    //
    unsigned long _neededFrames;

} LoudnessMeasurerChannel;


struct LoudnessMeasurer {
    size_t _channelCount;
    LoudnessMeasurerChannel *_channels;

    // How many samples fit in 100ms (rounded)
    unsigned long _samplesIn100ms;

    // How many samples fit in 400ms (rounded)
    unsigned long _samplesIn400ms;

    double _coefficients[2][5];
};


static inline void LoudnessMeasurerSetupFilter(LoudnessMeasurer *m, double sampleRate)
{
    double f0 = 1681.974450955533;
    double G  =    3.999843853973347;
    double Q  =    0.7071752369554196;

    double K  = tan(M_PI * f0 / sampleRate);
    double Vh = pow(10.0, G / 20.0);
    double Vb = pow(Vh, 0.4996667741545416);

    double pb[3] = {0.0,  0.0, 0.0};
    double pa[3] = {1.0,  0.0, 0.0};
    double rb[3] = {1.0, -2.0, 1.0};
    double ra[3] = {1.0,  0.0, 0.0};

    double a0 =      1.0 + K / Q + K * K      ;
    pb[0] =     (Vh + Vb * K / Q + K * K) / a0;
    pb[1] =           2.0 * (K * K -  Vh) / a0;
    pb[2] =     (Vh - Vb * K / Q + K * K) / a0;
    pa[1] =           2.0 * (K * K - 1.0) / a0;
    pa[2] =         (1.0 - K / Q + K * K) / a0;

    f0 = 38.13547087602444;
    Q  =  0.5003270373238773;
    K  = tan(M_PI * f0 / sampleRate);

    ra[1] =   2.0 * (K * K - 1.0) / (1.0 + K / Q + K * K);
    ra[2] = (1.0 - K / Q + K * K) / (1.0 + K / Q + K * K);

    m->_coefficients[0][0] = pb[0];
    m->_coefficients[0][1] = pb[1];
    m->_coefficients[0][2] = pb[2];
    m->_coefficients[0][3] = pa[1];
    m->_coefficients[0][4] = pa[2];

    m->_coefficients[1][0] = rb[0];
    m->_coefficients[1][1] = rb[1];
    m->_coefficients[1][2] = rb[2];
    m->_coefficients[1][3] = ra[1];
    m->_coefficients[1][4] = ra[2];
}


LoudnessMeasurer *LoudnessMeasurerCreate(unsigned int channelCount, double sampleRate, size_t totalFrames)
{
    LoudnessMeasurer *self = (LoudnessMeasurer *)calloc(1, sizeof(LoudnessMeasurer));
    if (!self) return NULL;

    self->_channelCount = channelCount;
    self->_channels = calloc(channelCount, sizeof(LoudnessMeasurerChannel));
    
    self->_samplesIn100ms = (sampleRate + 5) / 10;
    self->_samplesIn400ms = self->_samplesIn100ms * 4;
    
    for (int i = 0; i < channelCount; i++) {
        LoudnessMeasurerChannel *channel = &self->_channels[i];
        channel->_bufferFrames = self->_samplesIn400ms;

        channel->_blocksCapacity   = ceil(totalFrames / sampleRate) * 10;
        channel->_overviewCapacity = ceil(totalFrames / sampleRate) * 100;

        channel->_bufferPre  = (float  *)malloc( channel->_bufferFrames      * sizeof(float));
        channel->_overview   = (float  *)malloc( channel->_overviewCapacity  * sizeof(float));

        channel->_bufferPost = (double *)malloc( channel->_bufferFrames      * sizeof(double));
        channel->_scratch    = (double *)malloc((channel->_bufferFrames + 2) * sizeof(double));
        channel->_scratch2   = (double *)malloc((channel->_bufferFrames + 2) * sizeof(double));
        channel->_blocks     = (double *)malloc( channel->_blocksCapacity    * sizeof(double));

        channel->_neededFrames = self->_samplesIn400ms;
    }

    LoudnessMeasurerSetupFilter(self, sampleRate);

    return self;
}


void LoudnessMeasurerFree(LoudnessMeasurer *self)
{
    for (int c = 0; c < self->_channelCount; c++) {
        LoudnessMeasurerChannel *channel = &self->_channels[c];

        free(channel->_bufferPre);
        free(channel->_overview);

        free(channel->_bufferPost);
        free(channel->_scratch);
        free(channel->_scratch2);
        free(channel->_blocks);
    }

    free(self->_channels);

    free(self);
}


static inline void sBiquad(double *inBuffer, double *outBuffer, const double *coef, double *state, size_t frames)
{
    double *i = inBuffer  - 2;
    double *o = outBuffer - 2;

    i[0] = state[0];
    i[1] = state[1];
    o[0] = state[2];
    o[1] = state[3];

    vDSP_deq22D(i, 1, coef, o, 1, frames);

    state[0] = i[frames];
    state[1] = i[frames + 1];
    state[2] = o[frames];
    state[3] = o[frames + 1];
}


static inline void sFilter(const LoudnessMeasurer *self, LoudnessMeasurerChannel *channel, const float* src, size_t stride, size_t frames)
{
    if (frames == 0) return;

    float max = 0;
    vDSP_maxv(src, stride, &max, frames);

    float min = 0;
    vDSP_minv(src, stride, &min, frames);

    if (-min > max) max = -min;

    if (max > channel->_samplePeak) {
        channel->_samplePeak = max;
    }

#if __x86__
    unsigned int mxcsr = _mm_getcsr();
    _mm_setcsr(mxcsr | _MM_FLUSH_ZERO_ON);
#endif

    double *s1 = channel->_scratch  + 2;
    double *s2 = channel->_scratch2 + 2;

    vDSP_vspdp(src, stride, s1, 1, frames);

    float *bufferPre = channel->_bufferPre + channel->_bufferIndex;

    for (size_t i = 0; i < frames; i++) {
        bufferPre[i] = s1[i];
    }

    sBiquad(s1, s2, self->_coefficients[0], channel->_filterState[0], frames);
    sBiquad(s2, s1, self->_coefficients[1], channel->_filterState[1], frames);

    double *bufferPost = channel->_bufferPost + channel->_bufferIndex;

    for (size_t i = 0; i < frames; i++) {
        bufferPost[i] = s1[i];
    }
#if __x86__
    _mm_setcsr(mxcsr);
#endif
}


static inline void sCalculateOverview(const LoudnessMeasurer *self, LoudnessMeasurerChannel *channel, size_t framesPerBlock, size_t blockCount)
{
    size_t b = channel->_overviewCount == 0 ? 0 : blockCount - 10;

    for ( ; b < blockCount; b++) {
        size_t bufferIndex = channel->_bufferIndex + (framesPerBlock * b);
        bufferIndex = bufferIndex % channel->_bufferFrames;

        float max = 0;
        float m;
        
        if (channel->_overviewCount == 0 && bufferIndex == 0) {
            continue;
        }
        
        if (bufferIndex < framesPerBlock) {
            size_t i = 0;
            size_t frames = bufferIndex;

            if (frames > 0) {
                vDSP_maxv(&channel->_bufferPre[i], 1, &m, frames);
                if ( m > max) max =  m;

                vDSP_minv(&channel->_bufferPre[i], 1, &m, frames);
                if (-m > max) max = -m;
            }
            
            i = channel->_bufferFrames - (framesPerBlock - bufferIndex);
            frames = channel->_bufferFrames - i;

            vDSP_maxv(&channel->_bufferPre[i], 1, &m, frames);
            if ( m > max) max =  m;

            vDSP_minv(&channel->_bufferPre[i], 1, &m, frames);
            if (-m > max) max = -m;

        } else {
            size_t i      = bufferIndex - framesPerBlock;
            size_t frames = bufferIndex - i;

            if (bufferIndex == 0) {
                i = bufferIndex;
                frames = framesPerBlock;
            }

            vDSP_maxv(&channel->_bufferPre[i], 1, &m, frames);
            if ( m > max) max =  m;

            vDSP_minv(&channel->_bufferPre[i], 1, &m, frames);
            if (-m > max) max = -m;
        }

        size_t count = channel->_overviewCount;

        if (count < channel->_overviewCapacity) {
            channel->_overview[count] = max;
            channel->_overviewCount++;
        }
    }
}


static inline void sCalculateGatingBlock(const LoudnessMeasurer *self, LoudnessMeasurerChannel *channel, size_t framesPerBlock)
{
    double channelSum = 0.0;

    if (channel->_bufferIndex < framesPerBlock) {
        size_t i = 0;
        size_t frames =  channel->_bufferIndex;
        double acc = 0;
    
        if (frames) {
            vDSP_vsqD(&channel->_bufferPost[i], 1, channel->_scratch, 1, frames);
            vDSP_sveD(channel->_scratch, 1, &acc, frames);
        }
        channelSum += acc;

        i = channel->_bufferFrames - (framesPerBlock - channel->_bufferIndex);
        frames = channel->_bufferFrames - i;
        acc = 0;
        
        if (frames) {
            vDSP_vsqD(&channel->_bufferPost[i], 1, channel->_scratch, 1, frames);
            vDSP_sveD(channel->_scratch, 1, &acc, frames);
            channelSum += acc;
        }

    } else {
        size_t i      = channel->_bufferIndex - framesPerBlock;
        size_t frames = channel->_bufferIndex - i;

        if (frames) {
            vDSP_vsqD(&channel->_bufferPost[i], 1, channel->_scratch, 1, frames);
            vDSP_sveD(channel->_scratch, 1, &channelSum, frames);
        }
    }

    size_t count = channel->_blocksCount;

    if (count < channel->_blocksCapacity) {
        channel->_blocks[count] = channelSum;
        channel->_blocksCount++;
    }
}


static void sProcess(const LoudnessMeasurer *self, LoudnessMeasurerChannel *channel, const float *source, size_t stride, size_t inFrames)
{
    size_t index = 0;
    size_t frames = inFrames;

    while (frames > 0) {
        if (frames >= channel->_neededFrames) {
            sFilter(self, channel, source + (index * stride), stride, channel->_neededFrames);
            index  += channel->_neededFrames;
            frames -= channel->_neededFrames;
            channel->_bufferIndex += channel->_neededFrames;
            
            sCalculateGatingBlock(self, channel, self->_samplesIn400ms);
            sCalculateOverview(self, channel, self->_samplesIn100ms / 10, 40);

            // 100ms are needed for all blocks besides the first one
            channel->_neededFrames = self->_samplesIn100ms;

            if (channel->_bufferIndex == channel->_bufferFrames) {
                channel->_bufferIndex = 0;
            }

        } else {
            sFilter(self, channel, source + (index * stride), stride, frames);

            channel->_bufferIndex  += frames;
            channel->_neededFrames -= frames;

            frames = 0;
        }
    }
}


void LoudnessMeasurerScanAudioBuffer(LoudnessMeasurer *self, AudioBufferList *bufferList, size_t inFrames)
{
    dispatch_apply(self->_channelCount, dispatch_get_global_queue(0, 0), ^(size_t c) {
        LoudnessMeasurerChannel *channel = &self->_channels[c];
        sProcess(self, channel, bufferList->mBuffers[c].mData, 1, inFrames);
    });
}


NSData *LoudnessMeasurerGetOverview(LoudnessMeasurer *self)
{
    size_t  overviewCount = self->_channels[0]._overviewCount;
    UInt8  *overview      = malloc(sizeof(UInt8) * overviewCount);
    
    for (int i = 0; i < overviewCount; i++) {
        float max = 0;

        for (int c = 0; c < self->_channelCount; c++) {
            float m = self->_channels[c]._overview[i];
            if (m > max) max = m;
        }
    
        SInt16 result = floor(max * 255.0);
        
        if (result > 255) result = 255;
        if (result < 0)   result = 0;

        overview[i] = result;
    }
    
    return [[NSData alloc] initWithBytesNoCopy:overview length:(overviewCount * sizeof(UInt8)) freeWhenDone:YES];
}


double LoudnessMeasurerGetLoudness(LoudnessMeasurer *self)
{
    static double sRelativeGateFactor = 0;
    static double sAbsoluteGateFactor = 0;

    if (!sRelativeGateFactor) sRelativeGateFactor = pow(10.0, -10 / 10.0);
    if (!sAbsoluteGateFactor) sAbsoluteGateFactor = pow(10.0, (-70.0 + 0.691) / 10.0);

    double relativeThreshold = 0.0;
    double gatedLoudness = 0.0;
    size_t aboveThresholdCounter = 0;

    size_t  blocksCapacity = self->_channels[0]._blocksCapacity;
    size_t  blocksCount    = 0;
    double *blocks         = malloc(sizeof(double) * blocksCapacity);
    
    for (int i = 0; i < self->_channels[0]._blocksCount; i++) {
        double sum = 0;

        for (int c = 0; c < self->_channelCount; c++) {
            sum += self->_channels[c]._blocks[i];
        }

        sum /= (double)self->_samplesIn400ms;

        if (sum >= sAbsoluteGateFactor) {
            if (blocksCount < blocksCapacity) {
                blocks[blocksCount++] = sum;
            }
        }
    }


    for (int i = 0; i < blocksCount; i++) {
        aboveThresholdCounter++;
        relativeThreshold += blocks[i];
    }


    if (!aboveThresholdCounter) {
        free(blocks);
        return 0;
    }

    relativeThreshold /= (double) aboveThresholdCounter;
    relativeThreshold *= sRelativeGateFactor;
    aboveThresholdCounter = 0;

    for (int i = 0; i < blocksCount; i++) {
        double sum = blocks[i];

        if (sum >= relativeThreshold) {
            ++aboveThresholdCounter;
            gatedLoudness += sum;
        }
    }

    free(blocks);

    if (!aboveThresholdCounter) {
        return 0;
    }

    gatedLoudness /= (double)aboveThresholdCounter;
    return 10 * (log(gatedLoudness) / log(10.0)) - 0.691;
}


double LoudnessMeasurerGetPeak(LoudnessMeasurer *self)
{
    double result = 0;

    for (int c = 0; c < self->_channelCount; c++) {
        double peak = self->_channels[c]._samplePeak;
        if (peak > result) {
            result = peak;
        }
    }

    return result;
}

