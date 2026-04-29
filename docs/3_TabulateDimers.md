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
### Dependencies
Python modules pandas v1.4 and numpy v1.24

### Input
Text files output by MFEprimer dimer.

### Python syntax
```
import sys
sys.path.append('/multiplex_wormhole/src/')
from scripts.tabulate_MFEprimer_dimers import main as tabulateDimers
tabulateDimers(ALL_DIMERS, END_DIMERS, OUTPATH, OUTPRIMERPATH="False", deltaG=False)
```

### Command line syntax
```
cd multiplex_wormhole/src/scripts/
python3 tabulate_MFEprimer_dimers.py -a ALL_DIMERS.txt -e END_DIMERS.txt -o OUTPATH [-p PRIMERS_OUTPATH] [-d]
```

### Arguments
**ALL_DIMERS (-a)** : Text file output for all dimers output from MFE primer.

**END_DIMERS (-e)** : Text file output for end dimers output from MFE primer. If you did not predict end dimers (not recommended), you can provide a blank text file here.

**OUTPATH (-o)** : Directory path and prefix for output tables of primer pair interactions.

**OUTPRIMERPATH (-p)** : Directory path and prefix for output tables of primer interactions. (Default: "False" -produces no output)

**deltaG (-d)** : Tabulate dimer load using deltaG (True) or standard tally (False) method. (Default: False)

## Outputs
1. A N_PRIMER_PAIRS x N_PRIMER_PAIRS table which includes the dimer load estimated for all pairwise combinations of primer pairs.
2. A N_PRIMERS_PAIRS X N_PRIMER_PAIRS table which specifies pairwise primer pair interactions in binary (0 = no secondary structures between primer 1 and primer 2, 1 = at least 1 secondary structure between primer 1 and primer 2).
3. A N_PRIMERS x 2 table which includes the total number of primer dimers contributed by each primer pair.
4. A N_PRIMERS X 2 table which inclues the total number of primer pairs with which each primer pair interacts.

Optionally, tables can also be output for per-primer interactions between individual primers (not just pairs) if the `OUT_PRIMERPATH` argument is provided.

### Recommended usage in next steps
If you want to prioritize minimizing the number of pairwise interactions in a multiplex, use the binary outputs in the optimziation step. 

If you want to prioritize minimizing the total dimer load in a multiplex, use the count outputs in the optimization step.


[Previous](2_DimerPrediction.md)		[Next](4A_ExploreOptimParameters.md)
