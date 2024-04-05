# Filtering Primer Pairs with `filter_primers`

## Purpose
Primer pairs are filtered to avoid within-pair secondary structures (i.e., hairpins, homodimers, and heterodimers) based on primer3 output. 

By default, primer pairs are discarded if the following criteria are met:
- secondary structure Tm > 45 C

**AND**

- delta G < *threshold* where the threshold varies depending on the structure type:
  - hairpins: -2 kcal/mol
  - homodimers or heterodimers at primer ends: -5 kcal/mol
  - homodimers or heterodimers not at primer ends: -10 kcal/mol
 

## Usage

### Python syntax
`import os`

`os.chdir('/multiplex_wormhole')`

`from scripts.filter_primers import main as filterPrimers`

`filterPrimers(PRIMER_DIR, OUTPATH)`

with optional parameters:

`filterPrimers(PRIMER_DIR, OUTPATH, Tm_LIMIT=45, dG_HAIRPINS=-2000, dG_END_LIMIT=-5000, dG_MID_LIMIT=-10000)`


### Command line syntax
`cd multiplex_wormhole/scripts`

`python3 filter_primers.py PRIMER_DIR OUTPATH Tm_LIMIT dG_HAIRPINS dG_END_LIMIT dG_MID_LIMIT`


### Arguments
**PRIMER_DIR** : Path to directory holding primer3 .out output files created during batch primer design.

**OUTPATH** : Path where outputs will be saved (directory path + file prefix)

**Tm_LIMIT** (Optional): Maximum melting temperature allowed for secondary structures. Primer pairs with secondary structures that have melting temperatures above this value will be discarded. (Default: 45 Celsius)

**dG_HAIRPINS** (Optional): Minimum Gibbs free energy value (delta G) allowed for hairpin structures. Primers forming hairpins with delta G values lower than this value will be discarded. (Default: -2000 cal/mol = -2 kcal/mol)

**dG_END_LIMIT** (Optional): Minimum delta G allowed for primer dimers with binding on the primer ends. Primers forming end dimers with delta G lower than this value will be discarded. (Default: -5000 cal/mol = -5 kcal/mol)

**dG_MID_LIMIT** (Optional): Minimum delta G allowed for all other primer dimers. Primers forming dimers with delta G lower than this value will be discarded. (Default: -10000 cal/mol = -10 kcal/mol)


## Outputs
Primers that pass filtering are output to a CSV file (`OUTPATH`.csv) that includes one row per primer with the following information from the primer 3 output: 1. primerID
2. corresponding template sequenceID
3. primer pair number (N)
4. primer direction (FW or REV)
5. primer positing on template sequence (start BP)
6. primer length
7. annealing temperature (Tm)
8. primer proportion bound
9. amplicon size. 

PrimerIDs are defined by templateID.pairN.FW/REV. For example, the primerID "CLocus_704.5.FW" corresponds to the forward primer for primer pair 5 for template locus CLocus_704.

A list of template sequence IDs that had >1 primer pair(s) pass filtering are output to a text file named `OUTPATH`_LocusIDs.txt
