#!/bin/sh

## Change $BRiQ_BINPATH and $BRiQ_DATAPATH to the directories containing the
## binary executables and precompiled data files correspondingly.
export BRiQ_BINPATH=$HOME/RNA-BRiQ/build/bin
export BRiQ_DATAPATH=$HOME/BRiQ_data

INPUT=input     # Input file
OUTPDB=pred.pdb # Output: refined structure
RANDOMSEED=123  # Random seed

# RNA structure refinement
$BRiQ_BINPATH/BRiQ_Refinement $INPUT $OUTPDB $RANDOMSEED
