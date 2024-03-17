# Primer Prediction with <MFEprimer dimer>

## Purpose
Primer dimers are predicted using MFEprimer, with one MFEprimer output including all primer dimers expected to form and a second output including only primer dimers forming on the 3' end of primers. The following defaults are used for calculating dimer formation:

## Usage


## Defaults
These values are specified in the <multiplex_primer_design> script and may be altered. Specifications include:
- delta G threshold for structures at the 3' end: -3 kcal/mol
- delta G threshold for any other structures: -6 kcal/mol
- score limit: 3 (score calculated with +1 for each bp match and -1 for each bp mismatch)
- max mismatches in dimer: 40 bp
- dNTP concentration: 0.25 mM
- oligo concentration: 50 nM
- divalent cation concentration: 3.8 mM
- monovalent cation concentration: 50 mM
If arguments are excluded, MFEprimer's own defaults will be used, which may not match multiplex_wormhole defaults.
