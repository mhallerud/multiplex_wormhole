# multiplex_wormhole
*In silico* optimization for multiplex PCR assays that minimized predicted primer dimer loads. The pipeline was developed for genotyping by amplicon sequencing (i.e., reduced SNP panel) applications, however, the process is transferable to any multiplex PCR targeted sequencing approach. The impetus for multiplex_wormhole was genotyping noninvasive wildlife genetic samples. Default primer design settings are therefore conservative and tailored towards amplifying low concentration and degraded DNA such as that found in fecal and hair samples. 

Documentation webpage: [https://github.io/mhallerud/multiplex_wormhole]( https://github.io/mhallerud/multiplex_wormhole)

## Installation

## Dependencies
- Primer3 is used for primer design and can be downloaded [here for Windows](https://github.com/primer3-org/primer3/releases), or [here for Mac/Linux users](https://github.com/primer3-org/primer3). Keep track of where this downloads! The path to primer3_core will need to be updated at line 55 in multiplex_primer_design.py and line 88 in multiplex_wormhole.py.
- MFEprimer is used for dimer calculations and can be downloaded [here](https://www.mfeprimer.com/mfeprimer-3.1/#2-command-line-version). The path to MFEprimer-*-awd will need to be updated on line 54 of multiplex_primer_design.py and line 87 in multiplex_wormhole.py. MFEprimer Version 3.2.7 was used during development. 

The following Python packages are required and can normally be installed by running, e.g., "pip install pandas" on the command line or terminal:
- pandas (v1.4.4 used in development)
- numpy (v1.24.4 used in development)
- matplotlib (v3.5.2 used in development)
General Python modules required (these come pre-installed with most Python installations): os, sys, csv, random, math, signal, gc, itertools, shutil, glob, datetime

multiplex_wormhole was built and tested on MacOS with Python v3.9.13 in the Spyder IDE managed under Anaconda-Navigator. For those new to Python, [Anaconda](https://www.anaconda.com/products/navigator) is the recommended virtual environment manager. 


## Example Workflow
See [multiplex_primer_design](multiplex_primer_design.py) to run the full pipeline step by step.

## Comments/Questions/Ideas
Please report any problems or potential enhancements in the GitHub Issues page.