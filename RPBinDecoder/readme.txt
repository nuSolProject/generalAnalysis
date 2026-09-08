RPBinDecoderV1.cc for binVolts.c RP script data. General waveform reconstruction with full/gates integral and 3-exp fit for 2 channel data. This is V1 and it's decently stable, 57Co DP half-life is like 2 ns off so lets spend a couple months to fix that

RPBinDecoderV0.cc is the new and unimproved version of RPBinDecoder, but I spent like 5 minutes deleting the fit stuff and like 30 minutes writing comments. It's sorta tested, but not to the extent of V1. Simplfied main by removing input file loop timestamp shows integer precision in root file while nominal in coutRPCHData helper. Check stod-root ntpule compatibility or conversion function

These are just root macros. run with "root RPBinDecoderV0.cc" and I usually run with "-l -b -q" options for large batch runs
