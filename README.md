LoudnessMeasurer
================

OS X/iOS implementation of the EBU R128 loudness standard

This is heavily based on the excellent [libebur128](https://github.com/jiixyj/libebur128) and refactored such that:

1. All math uses Accelerate.framework (which in turn uses CPU vectorization)
2. Each channel is computed on a seperate thread via libdispatch

Note that:

3. LoudnessMeasurer only works on 2-channel files.
4. Only gated loudness and peak are calculated (no support for true peak)

That said, it provides an example of how you could speed up libebur128 using Accelerate.
