# Basic sequence alignment
This project implements simple global and local alignment algorithms that use naive and affine gap penalties in **C**.

Makefile uses **Makefile** and **gcc** as compiler. These packages should be installed into the system in order to run the project.
After navigating to the projects directory to compile the project execute
```
make
```

### Scoring Matrix for match/mismatch score

||A|C|G|T|
|-|-|-|-|-|
|A|2|-3|-3|-3|
|C|-3|2|-3|-3|
|G|-3|-3|2|-3|
|T|-3|-3|-3|2|

### Input
An input file is required. An example *sequences.fasta* file can be found in *test* folder.

### Parameters
- **--mode:** It will be selected from one of the followings:
  - **global:** Needleman-Wunsch with naive gap scoring
  - **local:** Smith-Waterman with naive gap scoring
  - **aglobal:** Needleman-Wunsch with affine gap scoring
  - **alocal:** Smith-Waterman with affine gap scoring
- **--input:** Input FASTA file for sequences
- **--gapopen:** Gap opening penalty for affine gap model, or unit gap cost for naive model
- **--gapext:** Gap extension penalty for affine gap model

### Output
- **global-naiveGap.aln** will be the only output file if the parameter is **--mode global**
- **global-affineGap.aln** will be the only output file if the parameter is **--mode aglobal**
- **local-naiveGap.aln** will be the only output file if the parameter is **--mode local**
- **local-affineGap.aln** will be the only output file if the parameter is **--mode alocal**

### Command Line Examples
Each line represents different runs for different modes. Be watchful of the options **--mode**, **--input**, **--gapopen**, **--gapext**.
```
allalign --mode global --input sequences.fasta --gapopen -5
allalign --mode aglobal --input sequences.fasta --gapopen -5 --gapext -2
allalign --mode local --input sequences.fasta --gapopen -5
allalign --mode alocal --input sequences.fasta --gapopen -5 --gapext -2
```
