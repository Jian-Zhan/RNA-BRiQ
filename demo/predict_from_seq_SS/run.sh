#!/bin/sh

## Change $BRiQ_BINPATH and $BRiQ_DATAPATH to the directories containing the
## binary executables and precompiled data files correspondingly.
export BRiQ_BINPATH=$HOME/RNA-BRiQ/build/bin
export BRiQ_DATAPATH=$HOME/BRiQ_data

INPUT=input     # Input file
OUTPDB=pred.pdb # Output: refined structure
RANDOMSEED=123  # Random seed

## Generate an initial PDB structure from the given sequence
$BRiQ_BINPATH/BRiQ_InitPDB seq.txt init.pdb

## Prepare the input file
echo "pdb init.pdb" > $INPUT
echo -n "seq " >> $INPUT
cat seq.txt >> $INPUT
echo -n "sec " >> $INPUT
cat ss.txt >> $INPUT
echo -n "nwc " >> $INPUT
cat nwc.txt >> $INPUT

## Predict RNA structure
$BRiQ_BINPATH/BRiQ_Predict $INPUT $OUTPDB $RANDOMSEED
