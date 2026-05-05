---
title: Panel Assessment
layout: page
permalink: /assess-panel
---
# Panel Assessment
The panel assessment function calculates the dimer load (both dimer count and mean deltaG) of an existing set of multiplex primers. The function runs [Dimer Prediction](2_DimerPrediction.md), then [Tabulates Dimers](3_TabulateDimers.md) for deltaG=False (i.e., counts) and deltaG=True (i.e., mean deltaG), and prints out a summary.


## Usage
### Python syntax
```
import multiplex_wormhole as mw
mw.assessPanel(PRIMERS,  ALL_DIMERS_dG=-8, END_DIMERS_dG=-4, BAD_DIMERS_dG=-10)
```

### Command line syntax
```
cd ~/multiplex_wormhole #navigate to where mw functions live
python3 panel_assessment.py  -i INPUT [-a ALLDIMERS_DG] [-e ENDDIMERS_DG] [-b BADDIMERS_DG]
```

### Arguments
**INPUT (-i)** : Input primers, can be either CSV or FASTA format. If CSV, must have "PrimerID" and "Sequence" fields for conversion.
**ALL_DIMERS_dG (-a)** : Gibbs free energy (delta G) threshold used to consider dimers not at 3' ends. [Default: -8]
**END_DIMERS_dG (-e)** : Gibbs free energy (delta G) threshold used to consider dimers at 3' ends. [Default: -4]
**BAD_DIMERS_dG (-b)** : Gibbs free energy (delta G) threshold used to consider dimers especially "bad". [Default: -10]


## Outputs
All outputs take the basename of the input primer file.
* `PREFIX`_MFEdimers.txt
* `PREFIX`_MFEdimers_ends.txt
* `PREFIX`_PrimerPairDimers_ + sum.csv, wide.csv, binary_sum.csv, binary_wide.csv
* `PREFIX`_PrimerPairDeltaG_ + mean.csv, wide.csv

Prints out:
* # primer pairs
* total # pairwise dimers
* # primer pairs involved in dimers
* dimers per pair
* mean deltaG of dimers
* Number of pairwise interactions with "bad" dimers


[Previous](4_OptimizeMultiplexPrimerSet.md)      [Next](8_SpecificityChecks.md)
