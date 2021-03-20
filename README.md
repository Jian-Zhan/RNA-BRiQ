# RNA-BRiQ

This repository contains codes for the following work: Peng Xiong, Ruibo Wu, Jian Zhan and Yaoqi Zhou, Robust RNA structure refinement by a nucleobase-centric sampling algorithm coupled with a backbone rotameric and quantum-mechanical-energy-scaled base-base knowledge-based potential.

This is free software for non-commercial users. If you have any questions, please contact Peng Xiong via email pengx_a@163.com


## Installation

### Prerequisites

Make sure these following prerequisites have been installed:
* C++ compiler: Needs to be C++11 compatible. Gcc 9.3.0 and Clang 9.1.0 tested.
* GNU Make
* Git
* CMake
* Wget and Tar - Optional: to download and extract the preprocessed data.

### Compilation

1. Clone from the repository:

```
git clone https://github.com/Jian-Zhan/RNA-BRiQ RNA-BRiQ
```

2. Change directory to `RNA-BRiQ/build/` and compile the codes:
```
mkdir RNA-BRiQ/build/
cd RNA-BRiQ/build/
cmake ../
make
```

### Data

Preprocessed data for the RNA-BRiQ can be downloaded from [http://servers.sparks-lab.org/downloads/RNA-BRiQ-data.tar.gz](http://servers.sparks-lab.org/downloads/RNA-BRiQ-data.tar.gz). The data should be downloaded, and extracted to `$BRiQ_DATAPATH` before continuing.

```
wget http://servers.sparks-lab.org/downloads/RNA-BRiQ-data.tar.gz

export BRiQ_DATAPATH=FILEPATH/BRiQ_data     ## Change "FILEPATH" to the real path to contain the preprocessed data
tar -xvz RNA-BRiQ-data.tar.gz --one-top-level=$BRiQ_DATAPATH
```

## Running the program 

### Environment variables

The path containing BRiQ executables should be appended to `$PATH`. And `$BRiQ_DATAPATH` should be correctly set.

```
export BRiQ_DATAPATH=FILEPATH/BRiQ_data    ## Change "FILEPATH" to the real path containing the preprocessed data
export PATH=$PATH:FILEPATH/BRiQ/build/bin  ## Change "FILEPATH" to the real path containing the compiled codes
```

### Input files

Before running the program, you need to prepare the input files. An example files can be found at `demo/gcaa/input`

```
pdb gcaa.pdb    ## Initial structure in PDB format
seq GCGCAAGC    ## RNA nucleotide sequence
sec ((....))    ## Watson-Crick pairs
nwc ..(..)..    ## Non-Watson-Crick pairs
fixed 0 1 6 7   ## Index of the fixed nucleotides during structure sampling, starting from 0
```

If there is only RNA sequence but no reference PDB structure, `BRiQ_init` can generate an initial PDB structure from sequence:

```
BRiQ_init GCGCAAGC init.pdb
```

And the base pairing information can be extracted by `BRiQ_assignSS` from the input PDB structure:

```
BRiQ_assignSS $INPUTPDB $OUTFILE
```

### RNA structure refinement
`BRiQ_Refinement` is the program to refine the initial structure of a RNA.

```
BRiQ_Refinement $INPUT $OUTPDB $RANDOMSEED
```

### RNA structure prediction
`BRiQ_Predict` is the program to predict the structure of a RNA.
```
BRiQ_Predict $INPUT $OUTPDB $RANDOMSEED
```
