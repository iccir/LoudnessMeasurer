//
//  main.m
//  loudness
//
//  Created by Ricci Adams on 2014-02-04.
//  Copyright (c) 2014 Ricci Adams. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <AudioToolbox/AudioToolbox.h>

#import "LoudnessMeasurer.h"

static double sProcessAudioFile(NSURL *url)
{
    double loudness = 0;
    double peak = 0;

    ExtAudioFileRef audioFile = NULL;
    OSStatus err = noErr;

    // Open file
    if (err == noErr) {
        err = ExtAudioFileOpenURL((__bridge CFURLRef)url, &audioFile);
        if (err) NSLog(@"ExtAudioFileOpenURL: %ld", (long)err);
    }


    AudioStreamBasicDescription fileFormat = {0};
    UInt32 fileFormatSize = sizeof(fileFormat);

    if (err == noErr) {
        err = ExtAudioFileGetProperty(audioFile, kExtAudioFileProperty_FileDataFormat, &fileFormatSize, &fileFormat);
    }


    AudioStreamBasicDescription clientFormat = {0};

    if (err == noErr) {
        UInt32 channels = fileFormat.mChannelsPerFrame;
        
        clientFormat.mSampleRate       = fileFormat.mSampleRate;
        clientFormat.mFormatID         = kAudioFormatLinearPCM;
        clientFormat.mFormatFlags      = kAudioFormatFlagIsFloat | kAudioFormatFlagIsNonInterleaved | kAudioFormatFlagsNativeEndian | kAudioFormatFlagIsPacked;
        clientFormat.mBytesPerPacket   = sizeof(float);
        clientFormat.mFramesPerPacket  = 1;
        clientFormat.mBytesPerFrame    = clientFormat.mFramesPerPacket * clientFormat.mBytesPerPacket;
        clientFormat.mChannelsPerFrame = channels;
        clientFormat.mBitsPerChannel   = sizeof(float) * 8;

        err = ExtAudioFileSetProperty(audioFile, kExtAudioFileProperty_ClientDataFormat, sizeof(clientFormat), &clientFormat);
    }
    
    SInt64 fileLengthFrames = 0;
    UInt32 fileLengthFramesSize = sizeof(fileLengthFrames);
   
    if (err == noErr) {
        err = ExtAudioFileGetProperty(audioFile, kExtAudioFileProperty_FileLengthFrames, &fileLengthFramesSize, &fileLengthFrames);
    }
    
    NSInteger bytesTotal = 0;

    if (err == noErr) {
        NSInteger framesRemaining = fileLengthFrames;
        NSInteger bytesRemaining = framesRemaining * clientFormat.mBytesPerFrame;
        NSInteger bytesRead = 0;

        bytesTotal = bytesRemaining;

        LoudnessMeasurer *measurer = LoudnessMeasurerCreate(clientFormat.mChannelsPerFrame, clientFormat.mSampleRate, framesRemaining);

        AudioBufferList *fillBufferList = alloca(sizeof(AudioBufferList) * clientFormat.mChannelsPerFrame);
        fillBufferList->mNumberBuffers = clientFormat.mChannelsPerFrame;
        
        for (NSInteger i = 0; i < clientFormat.mChannelsPerFrame; i++) {
            fillBufferList->mBuffers[i].mNumberChannels = clientFormat.mChannelsPerFrame;
            fillBufferList->mBuffers[i].mDataByteSize = clientFormat.mBytesPerFrame * 4096 * 16;
            fillBufferList->mBuffers[i].mData = malloc(clientFormat.mBytesPerFrame  * 4096 * 16);
        }

        while (1 && (err == noErr)) {
            UInt32 frameCount = (UInt32)framesRemaining;
            err = ExtAudioFileRead(audioFile, &frameCount, fillBufferList);

            if (frameCount) {
                LoudnessMeasurerScanAudioBuffer(measurer, fillBufferList, frameCount);
            } else {
                break;
            }
            
            framesRemaining -= frameCount;
        
            bytesRead       += frameCount * clientFormat.mBytesPerFrame;
            bytesRemaining  -= frameCount * clientFormat.mBytesPerFrame;

            if (framesRemaining == 0) {
                break;
            }
        }

        for (NSInteger i = 0; i < clientFormat.mChannelsPerFrame; i++) {
            free(fillBufferList->mBuffers[i].mData);
        }
        
        loudness = LoudnessMeasurerGetLoudness(measurer);
        peak = LoudnessMeasurerGetPeak(measurer);

        LoudnessMeasurerFree(measurer);
    }

    if (audioFile) {
        ExtAudioFileDispose(audioFile);
    }

    return loudness;
}


static double gr[] = {
    -23.0,
    -33.0,
    -23.0,
    -23.0,
    -23.0,
    -23.0,
    -23.0,
    -23.0,
    -23.0
};

static double gre[] = {
    -2.2953556442089987e+01,
    -3.2959860397340044e+01,
    -2.2995899818255047e+01,
    -2.3035918615414182e+01,
    -2.2949997446096436e+01,
    -2.3017157781104373e+01,
    -2.3017157781104373e+01,
    -2.2980242495081757e+01,
    -2.3009077718930545e+01
};


int main(int argc, const char * argv[])
{
@autoreleasepool {
    fprintf(stderr, "Note: the tests do not have to pass with EXACT_PASSED.\n"
        "Passing these tests does not mean that LoudnessMeasurer is "
        "100%% EBU R 128 compliant!\n\n");

    NSString *testsPath = [NSString stringWithCString:argv[1] encoding:NSUTF8StringEncoding];
    testsPath = [testsPath stringByAppendingPathComponent:@"tests"];

    void (^runTest)(const char *, NSInteger) = ^(const char *filename, NSInteger i) {
        NSString *filenameString = [NSString stringWithCString:filename encoding:NSUTF8StringEncoding];
        NSString *filePath = [testsPath stringByAppendingPathComponent:filenameString];
        NSURL *url = [NSURL fileURLWithPath:filePath];
        
        if (![[NSFileManager defaultManager] fileExistsAtPath:filePath]) {
            NSLog(@"%s not found.  Download https://tech.ebu.ch/docs/testmaterial/ebu-loudness-test-setv03.zip and extract into the tests directory", filename);
            exit(1);
        }
        
        double actual = sProcessAudioFile(url);
        
        if (actual == actual) {
            char *passOrFail = (actual <= gr[i] + 0.1 && actual >= gr[i] - 0.1) ? "PASSED" : "FAILED";
            char *exact = (actual == gre[i]) ?  "EXACT_PASSED" : "EXACT_FAILED";

            printf("%s, %s - %s: %1.16e\n", passOrFail, exact, filename, actual);
        }
    };

    // We only do gated loudness and (non-true) peak on stereo files, so this is a subset
    // of the libebur128 test suite
    //
    // Tests will be within 0.1dB, but won't be exact due to floating point errors
    //
    runTest("seq-3341-1-16bit.wav", 0);
    runTest("seq-3341-2-16bit.wav", 1);
    runTest("seq-3341-3-16bit-v02.wav", 2);
    runTest("seq-3341-4-16bit-v02.wav", 3);
    runTest("seq-3341-5-16bit-v02.wav", 4);
    runTest("seq-3341-7_seq-3342-5-24bit.wav", 7);
    runTest("seq-3341-2011-8_seq-3342-6-24bit-v02.wav", 8);
}
    return 0;
}

