---
title: Add Keeplist to FASTA
layout: page
permalink: /add-keeplist-to-fasta
---
# Add Keeplist to FASTA
Merges FASTA of newly designed primers with keeplist primers, accounting for overlapping templates (as long as sequence names are equivalent between the two files).

## Usage
### Python usage
```
from add_keeplist_to_fasta import main as AddKeeplist2FASTA
AddKeeplist2FASTA(MAIN_FA, KEEPLIST_FA, OUTPATH=None)
```

### Command line usage
```
python3 add_keeplist_to_fasta.py -i MAIN_FA -k KEEPLIST_FA [-o OUTPATH]
```

### Arguments
**MAIN_FA (-i)** : FASTA of new primers that keeplist will be merged into. Any overlapping sequence names will be removed from this file.
**KEEPLIST_FA (-k)** : FASTA of keeplist primers.
**OUTPATH (-o)** : Filepath to output FASTA. [Default: MAIN_FA prefix + "_plusKeeplist.fa"]
