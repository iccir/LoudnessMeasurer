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

static void sProcessAudioFile(NSString *path, BOOL directoryMode)
{
    NSURL *url = [NSURL fileURLWithPath:path];

    double loudness = 0;
    double peak = 0;

    ExtAudioFileRef audioFile = NULL;
    OSStatus err = noErr;

    // Open file
    if (err == noErr) {
        err = ExtAudioFileOpenURL((__bridge CFURLRef)url, &audioFile);
        if (err && directoryMode) return;
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

    
    fprintf(stdout, "%s: %lf\n", [[path lastPathComponent] UTF8String], loudness);
}



int main(int argc, const char * argv[])
{
@autoreleasepool {
    if (argc != 2) {
        fprintf(stderr, "Usage: loudness <dir_or_file>\n");
        exit(1);
    }

    NSString *arg = [NSString stringWithCString:argv[1] encoding:NSUTF8StringEncoding];
    BOOL directory = NO;

    if (![[NSFileManager defaultManager] fileExistsAtPath:arg isDirectory:&directory]) {
        fprintf(stderr, "Error: no audio file found at %s\n", argv[1]);
        exit(2);
    }
    
    if (directory) {
        NSError *error = nil;
        NSArray *contents = [[NSFileManager defaultManager] contentsOfDirectoryAtPath:arg error:&error];

        if (error) {
            NSLog(@"%@", error);
            exit(2);
        }

        for (NSString *subpath in contents) {
            sProcessAudioFile([arg stringByAppendingPathComponent:subpath], YES);
        }

    } else {
        sProcessAudioFile(arg, NO);
    }

}
    return 0;
}


