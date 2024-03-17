# Batch Primer Design with `primer3_batch_design`

## Purpose
Primers are designed for each template sequence using primer3, including calculating secondary structures using Illumina i5 and i7 adapters as overhangs.

## Usage
For primer3_batch_design to run, `primer3.sh` must be found in the same folder and settings files for Primer3 must be in the multiplex_wormhole/settings folder under `Primer3_Base_NoSecondaryFilters.txt` and `Primer3_Broad_NoSecondaryFilters.txt`.

### Python syntax
`import os`
`os.chdir('/multiplex_wormhole')`
`from scripts.primer3_batch_design import main as primer3BatchDesign`
`primer3BatchDesign(IN_CSV, OUTDIR, PRIMER3_PATH)`

### Command line syntax
`cd multiplex_wormhole/scripts`
`primer3_batch_design IN_CSV OUTDIR PRIMER3_PATH`

### Arguments
`IN_CSV` : Path to CSV file containing DNA template sequences in the following format (including headers):
| SEQUENCE_ID | SEQUENCE_TEMPLATE | SEQUENCE_TARGET |
| CLocus_704 | TCAGAGAC... | 53,36 |
| ... | ... | ... |

`OUTDIR`
`PRIMER3_PATH`


## Defaults:
Narrow primer3 settings (see --- file):
- Annealing Temp: 52 C
- Amplicon size: 70-120 bp
- Primer size: 18-26 bp (optimal: 20)
- Primer Tm: 57-63 Celsius (optimal: 60)
- Primer GC content: 30-70% (optimal: 50)
- GC Clamp: 1
- Max End GC: 4
- Max Poly X: 4
- dNTP concentration: 0.25 mM
- template concentration: 50 nM
- divalent cation concentration: 3.8 mM
- monovalent cation concentration: 50 mM

If primers can't be found for the above settings, constraints are relaxed (see --- file):
    - Primer GC content: 20-80%
    - GC Clamp: 0
    - Max End GC: 5
    - Max Poly X: 5
