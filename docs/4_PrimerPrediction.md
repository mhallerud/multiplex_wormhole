# Primer Prediction with `MFEprimer dimer`

## Purpose
Primer dimers are predicted using MFEprimer. MFEprimer dimer is run twice: first for all primer dimers expected to form, and then for primer dimers forming specifically on the 3' end of primers. This allows different filtering parameters to be used for end dimers and other dimers during dimer prediction.

## Usage
### Python syntax
`import os`

`os.system(MFEprimer_PATH+" dimer -i "+INPUT+" -o "+ALL_DIMERS+" -d -8 -s 3 -m 50 --diva 3.8 --mono 50 --dntp 0.25 --oligo 50")`

`os.system(MFEprimer_PATH+" dimer -i "+INPUT+" -o "+END_DIMERS+" -d -5 -s 3 -m 70 --diva 3.8 --mono 50 --dntp 0.25 --oligo 50 -p")`

### Command line syntax
`<MFEprimer_PATH> dimer -i <INPUT> -o <ALL_DIMERS> -d -8 -s 3 -m 50 --diva 3.8 --mono 50 --dntp 0.25 --oligo 50`

`<MFEprimer_PATH> dimer -i <INPUT> -o <END_DIMERS> -d -5 -s 3 -m 70 --diva 3.8 --mono 50 --dntp 0.25 --oligo 50 -p`

### Arguments
**MFEprimer_PATH** : Filepath to MFEprimer run file (e.g., mfeprimer-3.2.7-darwin-10.6-amd64).

**INPUT** : Filepath to CSV file output from filter_primers or check_primer_specificity containing primer information.

**ALL_DIMERS** : Filepath for MFEprimer output containing information on predicted general dimers.

**END_DIMERS** : Filepath for MFEprimer output containing information on predicted end dimers.


## Outputs
This step has a text file output for each run of MFEprimer dimer. The output is saved to the outpath provided (ALL_DIMERS or END_DIMERS).


## MFEprimer parameters & Defaults used in multiplex_wormhole
These values are specified in the <multiplex_primer_design> script and may be altered. Specifications include:
`-d` : maximum delta G threshold for structures in kcal/mol (Default: all dimers: -6, end dimers: -3)

`-s` : score limit for including predicted dimers, calculated with +1 for each bp match and -1 for each bp mismatch (Default: 3)

`m` : max number of mismatches allowed in dimer (Default: all dimers: 50 bp, end dimers: 70 bp)

`--dntp` : DNA template concentration in reaction (Default: 0.25 mM)

`--oligo` : primer concentration in reaction (Default: 50 nM)

`--diva` : divalent cation concentration (Default: 3.8 mM)

`--mono` : monovalent cation concentration (Default: 50 mM)

`p` : only output primers with 3' end bind

If arguments are excluded, MFEprimer's own defaults will be used, which may not match multiplex_wormhole defaults.

See the [MFEprimer website](https://www.mfeprimer.com) for more details on MFE primer.


## Citation
Wang, K. et al., MFEprimer-3.0: quality control for PCR primers. Nucleic acids research 47, W610–W613 (2019).
