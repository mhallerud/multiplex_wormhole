---
title: Convert CSV to FASTA
layout: page
permalink: /csv-to-fasta
---
# Convert CSV to FASTA
Converts CSV containing primers (one per row) to a FASTA file. If the output is to be used in other multiplex wormhole functions, primer sequences must be in the format <sequenceID>.<#>.FWD & <sequenceID>.<#>.REV - e.g., MACA01.0.FWD and MACA01.0.REV signify the "0th" set of primers that amplify the "MACA01" template.


## Usage
### Python usage
```
import multiplex_wormhole as mw
mw.CSVtoFASTA(IN_CSV, OUT_FA, ID_FIELD="PrimerID", SEQ_FIELD="Sequence", ENCODING=sys.getfilesystemencoding())
```

### Command line usage
```
cd ~/multiplex_wormhole #navigate to where your mw scripts live
python3 CSVtoFASTA.py -i INCSV -o OUTFA [-p PRIMERIDFIELD] [-s SEQFIELD] [-e CSV_ENCODING]
```

### Arguments
**IN_CSV (-i)** : CSV containing primer sequences. Must have ID_FIELD and SEQ_FIELD.

**OUT_FA (-o)** : Filepath to output FASTA.

**ID_FIELD (-p)** : Column name for primer names. [Default: PrimerID]

**SEQ_FIELD (-s)** : Column name for primer sequences. [Default: Sequence]

**ENCODING (-e)** : Encoding of CSV file. [Default: system encoding]
