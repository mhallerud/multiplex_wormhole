---
layout: default
nav_order: 1
parent: index
permalink: /primer-design
---

# Batch Primer Design
Designs primers for each template sequence provided using primer3 (Untergasser et al. 2012). Filters for secondary structure formation within pairs based on primer3 output, where primer pairs are removed if the melting temperature AND deltaG thresholds are not met for hairpins, homodimers, and heterodimers. Outputs primer sequences and details to a table and FASTA file(s)- at this stage, keeplist primers are added to the FASTA.

**Important: By default, primers are designed with the Illumina Nextera i5 and i7 adapters attached. If other adapter sequences are desired, use SETTINGS={'SEQUENCE_OVERHANG_LEFT': "your_FWD_adapter_sequence", 'SEQUENCE_OVERHANG_RIGHT': "your_REV_adapter_sequence"}**


## Usage
### Dependencies
Requires primer3-py Python module

### Python syntax
```
mw.primer3BatchDesign(TEMPLATES, OUTPATH, Tm_LIMIT=45.0, dG_HAIRPINS=-2, dG_END_LIMIT=-4, dG_MID_LIMIT=-8, KEEPLIST=None, ENABLE_BROAD=False, SETTINGS=None)
```

### Command line syntax
```
cd ~/multiplex_wormhole/src/multiplex_wormhole
python3 batch_primer3_design.py -t TEMPLATES -o OUT [-l TM_LIMIT] [-d HAIRPINS_DG] [-m MIDDIMERS_DG] [-e ENDDIMERS_DG] [-k KEEPLIST] [-s SETTINGS] [-b]
```

### Arguments
**TEMPLATES (-t)** : CSV filepath to templates in the following format (including these specific headers):

| SEQUENCE_ID   | SEQUENCE_TEMPLATE    | SEQUENCE_TARGET    |
| ------------- | -------------------- | ------------------ |
| CLocus_704    | TCAGAGAC...          | 53,1               |
| ...           | ...                  | ...                |

The sequence target is in <POSITION,LENGTH> format. In the example, the target for PCR is a 1-bp region (i.e., SNP or indel) at basepair 53. (Required)

**OUTPATH (-o)** : Prefix (including directory structure) for output files which will include a CSV and FASTA. (Required)

**Tm_LIMIT (-l)** : Upper limit for melting temperatures (in Celsius) of intra-pair secondary structures, including hairpins, homodimers, and heterodimers. (Default: 45.0)

**dG_HAIRPINS (-d)** : Lower limit for delta G (i.e., Gibbs free energy) of hairpin structures. (Default: -2)

**dG_MID_LIMIT (-m)** : Lower limit for delta G (i.e., Gibbs free energy) of primer dimers not occurring at 3' ends. (Default: -8)

**dG_END_LIMIT (-e)** : Lower limit for delta G (i.e., Gibbs free energy) of primer dimers at 3' ends. (Default: -4)

**KEEPLIST (-k)** : FASTA of keeplist primers, where names follow the format <SeqID>.<#>.FWD and <SeqID>.<#>.REV (e.g., MACA01.1.FWD, MACA01.1.REV). *IMPORTANT*: If the same target is in the keeplist and the templates file, ensure that SeqIDs name match- otherwise, multiple primer pairs may be designed for the same target. (Default: None)

**ENABLE_BROAD (-b)** : Try less conservative settings if primers could not be developed for a template under strict design settings? (Default: False)

**SETTINGS (-s)** : Primer3 settings, provided in dictionary format (e.g. `{'PRIMER_ANNEALING_TEMP': 52.0, 'SEQUENCE_OVERHANG_LEFT': "tcgtcggcagcgt..."}`). If running from the command line, use the format `-s '{"PRIMER_ANNEALING_TEMP": "52.0"}'`. Setting names and definitions can be found on the [primer3 website](https://htmlpreview.github.io/?https://github.com/primer3-org/primer3/blob/v2.6.1/src/primer3_manual.htm#globalTags). Default settings are in lines 72-187 of the [primer3_batch_design.py script](https://github.com/mhallerud/multiplex_wormhole/blob/main/src/scripts/primer3_batch_design.py). (Default: None)

**(-h)**: Help for usage.


## Outputs
* `OUTPATH`.csv : CSV of filtered primer pair details output by primer3, including primerIDs, sequences, amplicon size, start BP within template, primer length, melting temperature, % GC content, proportion bound at annealing temp, 3' end stability, penalty, and template mispriming.
* `OUTPATH`.fa : FASTA of filtered primer pairs.
* `OUTPATH`_plusKeeplist.fa : FASTA merging filtered primer pairs with KEEPLIST primer pairs (automatically removes any duplicate names & primers between the KEEPLIST and TEMPLATES).


## Default Settings
Default settings are largely based on Eriksson et al. (2020). Important ones include:

* 'SEQUENCE_OVERHANG_LEFT'="tcgtcggcagcgtcagatgtgtataagagacag" #Illumina Nextera i5 (5'-->3')
* 'SEQUENCE_OVERHANG_RIGHT'="gtctcgtgggctcggagatgtgtataagagacag" #Illumina Nextera i7 (5'-->3')
* 'PRIMER_ANNEALING_TEMP'=52 #primer annealing temp (Celsius)
* 'PRIMER_PRODUCT_SIZE_RANGE'=70-120 #allowed amplicon size (bp)
* 'PRIMER_PRODUCT_OPT_SIZE'=100 #optimal amplicon size (bp)
* 'PRIMER_OPT_SIZE'=20 #optimal primer size'
* PRIMER_MIN_SIZE'=18 #min primer size
* 'PRIMER_MAX_SIZE'=26 #max primer size
* 'PRIMER_MIN_TM'=57 #min primer melting temp
* 'PRIMER_MAX_TM'=63 #max primer melting temp
* 'PRIMER_OPT_TM'=60 #optimal primer melting temp
* 'PRIMER_MAX_GC'=70 #max primer GC content
* 'PRIMER_MIN_GC'=30 #min primer GC content
* 'PRIMER_OPT_GC_PERCENT'=50 #optimal primer GC content
* 'PRIMER_MAX_END_GC'=4 #max of for GC at 3' end
* 'PRIMER_MAX_POLY_X'=4 #max nucleotide repeats
* 'PRIMER_DNTP_CONCENTRATION'=0.25 #primer concentrations in PCR (nM)
* 'PRIMER_DNA_CONCENTRATION'=50 #template DNA concentration (nM)
* 'PRIMER_SALT_DIVALENT'=3.8 #divalent salt cation concentration (mM)
* 'PRIMER_SALT_MONOVALENT'=50 #monovalent salt cation conc (mM)


## Citations
Eriksson, CE, Ruprecht J, Levi T. 2020. More affordable and effective noninvasive SNP genotyping using high-throughput amplicon sequencing. Molecular Ecology Resources 20(4): [doi: 10.1111/1755-0998.13208](https://doi.org/10.1111/1755-0998.13208).

Untergasser, A, Cutcutache, I, Koressaar, T, Ye, J, Faircloth, BC, Remm, M, Rozen SG. 2012. Primer3--new capabilities and interfaces. Nucleic Acids Research: e115. [doi: 10.1093/nar/gks596](https://doi.org/10.1093/nar/gks596). 

[Next](2_DimerPrediction.md)
