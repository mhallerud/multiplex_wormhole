---
title: Primer Dimer Prediction
layout: page
permalink: /dimer-prediction
nav_order: 1
parent: index
---
# Dimer Prediction
Primer dimers are predicted using [MFEprimer v3](https://www.mfeprimer.com) (Wang et al. 2019). The `MFEprimer dimer` function is run twice: first for all primer dimers expected to form, and then for primer dimers forming on the 3' end of primers which are more problematic. This allows different filtering parameters to be used for end dimers.


## Input
By default, multiplex_wormhole uses the FASTA output from the batch primer design step, with keeplist primers added. 

    **IMPORTANT**: All candidate primers need to be included at this step, including keeplist primers. These should have been added in the previous step, otherwise you can use add_keeplist_to_fasta.py. For highest accuracy, primers should include adapter sequences. Primers should be in the 5'-->3' orientation.


## Usage
### Dependencies
[MFEprimer3.2.7](https://github.com/quwubin/MFEprimer-3.0/releases/) must be downloaded and located in the /multiplex_wormhole/src directory. MFEprimer can be auto-downloaded with the setup_mfeprimer.py script- [see the homepage Quick Start](index.md#quick-start) for instructions, or otherwise manually downloaded and copied to /multiplex_wormhole/src.

### Python syntax
```
# setup
import sys
import glob
import subprocess
sys.path.append("/multiplex_wormhole/src")
MFEprimer_PATH = glob.glob("/multiplex_wormhole/src/*mfeprimer*)[0]

subprocess.call(MFEprimer_PATH+" dimer -i "+INPUT+" -o "+ALL_DIMERS+" -d -8 -s 3 -m 50 --diva 3.8 --mono 50 --dntp 0.25 --oligo 50", shell=True)
subprocess.call(MFEprimer_PATH+" dimer -i "+INPUT+" -o "+END_DIMERS+" -d -3 -s 3 -m 70 --diva 3.8 --mono 50 --dntp 0.25 --oligo 50 -p", shell=True) 
```

### Command line syntax
```
/multiplex_wormhole/src/*mfeprimer* dimer -i <INPUT> -o <ALL_DIMERS> -d -8 -s 3 -m 50 --diva 3.8 --mono 50 --dntp 0.25 --oligo 50
/multiplex_wormhole/src/*mfeprimer* dimer -i <INPUT> -o <END_DIMERS> -d -5 -s 3 -m 70 --diva 3.8 --mono 50 --dntp 0.25 --oligo 50 -p
```

### Parameters
**INPUT (-i)** : FASTA of primer sequences (with adapters) used to predict dimers. (Required)

**ALL_DIMERS (-o)** : Filepath for MFEprimer text output containing information on predicted dimers.

**END_DIMERS (-o)** : Filepath for MFEprimer text output containing information on predicted 3' end dimers.

**(-d)** : Maximum (i.e., least negative) deltaG (kcal/mol) threshold of predicted dimers.

**(-s)** : Minimum score threshold of predicted dimers, where scores are calculated with +1 for each match and -1 for each mismatch. 

**(-m)** : Max allowed mismatches per dimer.

**(--diva)** : Concentration of divalent salt cations (mM).

**(--mono)** : Concentration of monovalent salt cations (mM).

**(--dntp)** : Concentration of dNTPs (mM).

**(--oligo)** : Concentration of primers (nM).

**(-p)** : Only output dimers with 3' end bind?


## Outputs
Two textfiles saved at ALL_DIMERS and END_DIMERS, each of which contain binding information for all predicted dimers, e.g.:

Dimer 164: 17249.3.REV x 19580.2.REV
Score: 6, Delta G = -5.82 kcal/mol

gtctcgtgggctcggagatgtgtataagagacagttctctggccttcctccag
                                               ::::::
                                               gaggtcaacaaccactaccggacagagaatatgtgtagaggctcgggtgctctg


## Alternative dimer calculation tools
The pipeline is built to use MFEprimer dimer to calculate dimer formation, however the optimization process will accept any input tables as long as the 2 input tables specify 1) pairwise dimer loads between primer pairs and 2) the total dimer load per primer pair, with primer pair IDs matching between the input templates and both tables. 

One alternative to MFEprimer dimer is Primer Suite's Primer Dimer tool, though we have found dimers to be overestimated using this tool. To use this tool to calculate dimers:
1. Copy the *_specificityCheck_passed.csv file into [primer-dimer.com](https://primer-dimer.com), select 'Multiplex Analysis' and 'Dimer Structure Report', and click Submit. Depending on how many loci you provided, this step may take awhile (~20 minutes for primers for 50 loci, multiple hours for 1200 loci).

2. Run [translate_primerSuite_report.R](../scripts/extra/translate_primerSuite_report.R) on the PrimerDimerReport.txt file downloaded from primer-dimer.com to convert this to a CSV. The delta G threshold specified will filter out any primer dimers with delta G above this value.

Another alternative is [ThermoFisher's Multiplex Primer Design](https://www.thermofisher.com/us/en/home/brands/thermo-scientific/molecular-biology/molecular-biology-learning-center/molecular-biology-resource-library/thermo-scientific-web-tools/multiple-primer-analyzer.html), however multiplex_wormhole currently has no support for translating these outputs into tables.


## Citations
Wang, K. et al., MFEprimer-3.0: quality control for PCR primers. Nucleic acids research 47, W610–W613 (2019).


[Previous](1_BatchPrimerDesign.md)		[Next](3_TabulateDimers.md)
