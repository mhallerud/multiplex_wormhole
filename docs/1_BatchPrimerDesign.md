# Batch Primer Design with `primer3_batch_design`

## Purpose
Primers are designed for each template sequence using primer3, including predicting secondary structure formation for adapter-ligated primers.


## Usage
For primer3_batch_design to run, ***primer3.sh*** must be found in the same folder and settings files for Primer3 must be in the multiplex_wormhole/settings folder under ***Primer3_Base_NoSecondaryFilters.txt*** and ***Primer3_Broad_NoSecondaryFilters.txt***.

### Python syntax
`import os`

`os.chdir('/multiplex_wormhole')`

`from scripts.primer3_batch_design import main as primer3BatchDesign`

`primer3BatchDesign(IN_CSV, OUTDIR, PRIMER3_PATH)`

### Command line syntax
`cd multiplex_wormhole/scripts`

`python3 primer3_batch_design.py IN_CSV OUTDIR PRIMER3_PATH`

### Arguments
**IN_CSV** : Path to CSV file containing DNA template sequences in the following format (including headers):

| SEQUENCE_ID   | SEQUENCE_TEMPLATE    | SEQUENCE_TARGET    |
| ------------- | -------------------- | ------------------ |
| CLocus_704    | TCAGAGAC...          | 53,1               |
| ...           | ...                  | ...                |

The sequence target is in <POSITION,LENGTH> format. In the example, there is a SNP at basepair 53 within the locus 704.

**OUTDIR** : Directory where primer output files will be saved. This directory must already exist.

**PRIMER3_PATH** : Path to primer3_core file.



## Defaults
Default settings are found in the ***Primer3_Base_NoSecondaryFilters.txt*** and ***Primer3_Broad_NoSecondaryFilters.txt*** under multiplex_wormhole/primer3_settings. These settings can be manually changed in the text files. For details on primer3 setting options and definitions, see the [primer3 Manual](https://primer3.org/manual.html). Settings can also be explored in [Primer3Plus](https://www.primer3plus.com) and settings file saved. Make sure that SEQUENCE_ID=, SEQUENCE_TEMPLATE=, and SEQUENCE_TARGET= remain blank before inputting to the script.

These default settings follow Eriksson et al. 2020 and are intended for amplifying DNA for SNP-based genotyping assays of wildlife from degraded samples (specifically, noninvasive genetic samples such as scats). Specifically, some of the important defaults include:

Illumina Nextera i5 and i7 adapters are added to 5'ends of output primers
- SEQUENCE_OVERHANG_LEFT=tcgtcggcagcgtcagatgtgtataagagacag
- SEQUENCE_OVERHANG_RIGHT=gtctcgtgggctcggagatgtgtataagagacag

Primer annealing temp is 52 Celsius
- PRIMER_ANNEALING_TEMP=52

Amplicon size is 70-120 base pairs, with an optimal of 100 bp
- PRIMER_PRODUCT_SIZE_RANGE=70-120
- PRIMER_PRODUCT_OPT_SIZE=100

Primer size range is 18-26 bp (optimal 20 bp)
- PRIMER_OPT_SIZE=20
- PRIMER_MIN_SIZE=18
- PRIMER_MAX_SIZE=26

Primer Tm range from 57-63 (optimal 60)
- PRIMER_MIN_TM=57
- PRIMER_MAX_TM=63
- PRIMER_OPT_TM=60

Primer GC content must be between 30-70% (optimal 50%)
- PRIMER_MAX_GC=70
- PRIMER_MIN_GC=30
- PRIMER_OPT_GC_PERCENT=50

Primers must have at least one G or C (a GC clamp) at the 3' end
- PRIMER_GC_CLAMP=1

At most, primers can have 4 Gs or Gcs at the ends
- PRIMER_MAX_END_GC=4

At most, primers can have 4 repeats of the same base
- PRIMER_MAX_POLY_X=4

Primer concentrations are 0.25 nM in the final PCR reaction
- PRIMER_DNTP_CONCENTRATOIN=0.25

DNA template concentration is 50 nM in the final PCR reaction
- PRIMER_DNA_CONCENTRATION=50

Salt concentrations are 3.8 mM for divalent cations and 50 mM for monovalenet cations
- PRIMER_SALT_DIVALENT=3.8
- PRIMER_SALT_MONOVALENT=50

If primers can't be found for the above settings, constraints are relaxed:
- Primer GC content: 20-80%
- GC Clamp: 0
- Max End GC: 5
- Max Poly X: 5

## Outputs
For each template sequence, a Primer3 output file and error file will be saved with naming based on the template SEQUENCE_ID: <SEQUENCE_ID>.out and <SEQUENCE_ID>.err. All outputs will be saved to `OUTDIR`. See the [primer3 Manual](https://primer3.org/manual.html) for interpreting output files.


## Citations
Eriksson, CE, Ruprecht J, Levi T. 2020. More affordable and effective noninvasive SNP genotyping using high-throughput amplicon sequencing. Molecular Ecology Resources 20(4): 10.1111/1755-0998.13208.
