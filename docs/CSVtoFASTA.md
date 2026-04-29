---
title: Convert CSV to FASTA
layout: page
permalink: /assess-panel
---
# Convert CSV to FASTA
Converts CSV containing primers (one per row) to a FASTA file. If used in optimization functions, primer sequences must be in the format <sequence_name>.<#>.FWD & <sequence_name>.<#>.REV - e.g., MACA01.0.FWD and MACA01.0.REV signify the "0th" set of primers that amplify "MACA01" template.


## Usage
### Python usage
```
from scripts.CSVtoFasta import main as CSVtoFASTA
CSVtoFASTA(IN_CSV, OUT_FA, ID_FIELD="PrimerID", SEQ_FIELD="Sequence", ENCODING=sys.getfilesystemencoding())
```

### Command line usage
```
python3 CSVtoFASTA.py -i INCSV -o OUTFA [-p PRIMERIDFIELD] [-s SEQFIELD] [-e CSV_ENCODING]
```

### Arguments
**IN_CSV (-i)** : CSV containing primer sequences. Must have ID_FIELD and SEQ_FIELD.
**OUT_FA (-o)** : Filepath to output FASTA.
**ID_FIELD (-p)** : Column name for primer names. [Default: PrimerID]
**SEQ_FIELD (-s)** : Column name for primer sequences. [Default: Sequence]
**ENCODING (-e)** : Encoding of CSV file. [Default: system encoding]
