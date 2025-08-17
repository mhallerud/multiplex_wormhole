# Checking Primer Specificity with `check_primer_specificity`

## Purpose
Specificity of primers are checked against all templates. Any primer pairs that aren't specific to one sequence are discarded to avoid off-target amplification. 

When a genome for the target species is available, this step should be used in conjunction with pre-filtering target sequences to ensure that template sequences only match one location on the genome. This can be done via short-read alignment tools such as bwa mem or BLAST.

## Implementation

### Python syntax
```
import os
os.chdir('/multiplex_wormhole')
from scripts.check_primer_specificity import main as specificityCheck
specificityCheck(PRIMERS, TARGET, OUTPATH)
```

## Command line syntax
```
cd multiplex_wormhole/scripts
python3 check_primer_specificity.py PRIMERS TARGET OUTPATH
```

## Arguments
**PRIMERS** : Filepath of primer CSV that was output from the filtering step.

**TARGET** : Filepath to a templates CSV to check primer specificity against.

**OUTPATH** : Directory and file prefix for outputs.


## Outputs
**`OUTPATH`_passed.csv** : CSV file containing all primers passing the specificity check. File format is equivalent to the CSV output from `filter_primers`.

**`OUTPATH`_passed.fa** : Fasta file containing all primers passing the specificity check. Sequence IDs are primer IDs.

**`OUTPATH`_failed.csv** : CSV file containing primers that failed the specificity check. File format is equivalent to the other CSVs.


[Previous](2_FilterPrimers.md)		[Next](4_DimerPrediction.md)