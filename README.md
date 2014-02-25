LoudnessMeasurer
================

OS X/iOS implementation of the EBU R128 loudness standard

This is heavily based on the excellent (libebur128)[https://github.com/jiixyj/libebur128] but refactored such that:

1) All math uses Accelerate.framework (which in turn uses CPU vectorization)
2) Each channel is computed on a seperate thread via libdispatch
