---
title: Tabulate Dimers
layout: page
permalink: tabulate-dimers
nav_order: 1
parent: index
---
# Tabulate Dimers
Converts text files output by MFEprimer dimer into tables of pairwise primer pair interactions and total secondary interactions per primer pair. For pairwise interactions, dimer load can be summarized as the count of dimers between two pairs (standard mode) or as the minimum (i.e., worst) delta G between two primer pairs (deltaG mode). For each primer pair, interactivity can be summarized as the total count of dimers across all primer pairs (standard mode) or as the mean delta G across all combinations (deltaG mode).
   
## Usage
### Python syntax
```
mw.tabulateDimers(ALL_DIMERS, END_DIMERS, OUTPATH, OUTPRIMERPATH="False", deltaG=False)
```

### Command line syntax
```
cd ~/multiplex_wormhole/
python3 tabulate_dimers.py -a ALL_DIMERS.txt -e END_DIMERS.txt -o OUTPATH [-p PRIMERS_OUTPATH] [-d]
```

### Arguments
**ALL_DIMERS (-a)** : Text file output for all dimers output from MFE primer.

**END_DIMERS (-e)** : Text file output for end dimers output from MFE primer. If you did not predict end dimers (not recommended), you can provide a blank text file here.

**OUTPATH (-o)** : Directory path and prefix for output tables of primer pair interactions.

**OUTPRIMERPATH (-p)** : Directory path and prefix for output tables of primer interactions. (Default: "False" -produces no output)

**deltaG (-d)** : Tabulate dimer load using deltaG (True) or standard tally (False) method. (Default: False)

## Outputs
For deltaG=True:
* **`OUTPATH`_wide.csv** : N_PRIMER_PAIRS x N_PRIMER_PAIRS table with minimum deltaG per pairwise combination of primer pairs (i.e., the worst pairwise interaction).
* **`OUTPATH`_mean.csv** : N_PRIMER_PAIRS x 2 table with the mean deltaG per primer pair across all other primer pairs.

For deltaG=False:
* **`OUTPATH`_wide.csv** : N_PRIMER_PAIRS x N_PRIMER_PAIRS table with the count of dimers per pairwise combination of primer pairs.
* **`OUTPATH`_sum.csv** : N_PRIMER_PAIRS x 2 table with the sum of dimers for each predicted primer pair.
* **`OUTPATH`_binary_wide.csv** & **`OUTPATH`_binary_sum.csv** : Equivalent to the above, but with binary outputs for dimers rather than counts.

Equivalent outputs are provided for individual primers if the `OUTPRIMERPATH` field is provided.

### But which do I use moving forwards?
* If you want to prioritize minimizing the number of pairwise interactions in a multiplex, use the binary outputs in the optimization step. 
* If you want to prioritize minimizing the total dimer load in a multiplex, use the count outputs in the optimization step.
* If you want to prioritize minimizing the worst dimers in a multiplex (e.g., for complex problems where dimers are inevitable), use the deltaG outputs in the optimization step.


[Previous](2_DimerPrediction.md)		[Next](4A_ExploreOptimParameters.md)
