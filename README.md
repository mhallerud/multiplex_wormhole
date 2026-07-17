# multiplex_wormhole
*In silico* optimization for multiplex PCR assays that minimized predicted primer dimer loads. The pipeline was developed for genotyping by amplicon sequencing (i.e., reduced SNP panel) applications, however, the process is transferable to any multiplex PCR targeted sequencing approach. The impetus for multiplex_wormhole was genotyping noninvasive wildlife genetic samples. Default primer design settings are therefore conservative and tailored towards amplifying low concentration and degraded DNA such as that found in fecal and hair samples. 

Latest version on PyPi: [https://pypi.org/project/multiplex-wormhole/](https://pypi.org/project/multiplex-wormhole/)
Full documentation: [https://mhallerud.github.io/multiplex_wormhole](https://mhallerud.github.io/multiplex_wormhole/)


## Installation
### Set up a virtual environment
multiplex_wormhole was built and tested on MacOS with Python v3.9.13 in the Spyder IDE managed under Anaconda-Navigator. For those new to Python or with existing python packages, [Anaconda](https://www.anaconda.com/products/navigator) is the recommended virtual environment manager. For a conda virtual environment within your working directory:
```
conda create -n py39 python=3.9 #create new virtual env w/ python v3.9
conda activate py39 #enter virtual env
```

### Install multiplex wormhole
```
pip install multiplex_wormhole
multiplex-wormhole -h
```

### Or clone the repo
If pip install doesn't work, you can also install manually by taking the following steps (from the command line):
```
# clone GitHub repo
git clone https://github.com/mhallerud/multiplex_wormhole

# install dependencies
pip install -r multiplex_wormhole/requirements.txt

# install as python package
pip install -e multiplex_wormhole
```

Now you are ready to run multiplex wormhole!


## Quick Start
### Command line syntax

```
# PANEL DESIGN (minimum inputs)
# usage: multiplex-wormhole [-h] -t TEMPLATES -n NLOCI -o OUTDIR [-p PREFIX]
#                           [-k KEEPLIST] [-r RUNS] [-i ITER] [-c CYCLES]
#                           [-s GREEDY] [-d] [-v]
# example for standard optimization with minimum settings + prefix & keeplist:
multiplex-wormhole -t "Input_Templates.csv" -n 20 -o "Test_MW" -p "Test_MW_default" -k "Keeplist.fa"

# PANEL ASSESSMENT
# usage: mw-assess-panel [-h] -i INPUT [-a ALLDIMERS_DG] [-e ENDDIMERS_DG] [-b BADDIMERS_DG]
# example with defaults:
mw-assess-panel -i "Primers.fasta" -a -8 -e -4 -b -10
```

### Python syntax
```
# load module
import multiplex_wormhole as mw

# panel design example (showing defaults)
mw.multiplexWormhole(TEMPLATES="Input_Templates.csv", 
                     N_LOCI=50, 
                     OUTDIR="Test_MW", 
                     PREFIX="Test_MW_default",
                     KEEPLIST_FA="Keeplist.fa",
                     N_RUNS=10, ITERATIONS=1000, CYCLES=10, GREEDY=5000, deltaG=False, VERBOSE=False)#optional

# panel assessment example (showing defaults)
mw.assessPanel(PRIMERS,
               ALL_DIMERS_dG=-8, END_DIMERS_dG=-4, BAD_DIMERS_dG=-10) #optional
```

### Arguments

#### multiplex-wormhole

* **TEMPLATES (-t –-templates)** : Path to templates CSV. NLOCI (-n –nloci) : Final panel size (i.e., # primer pairs & # templates amplified). 

* **OUTDIR (-o –-outdir)** : Filepath where output directory will be created and all outputs saved within a generated folder structure. 

* **PREFIX (-p –-prefix)** : Prefix for all outputs. [Defaults to a timestamp if None provided] 

* **KEEPLIST_FA (-k –-keeplist)** : Path to keeplist FASTA. [Default: None] 

* **N_RUNS (-r –-runs)** : Number of optimization runs. [Default: 10] 

* **ITERATIONS (-i –-iter)** : Number of simulated annealing iterations per cycle. [Default: 1000] 

* **CYCLES (-c –-cycles)** : Number of simulated annealing cycles per run. [Default: 10] 

* **GREEDY (-s --simple)** : Number of simple iterative improvement iterations per run. [Default: 5000] 

* **deltaG (-d –-deltaG)** : Optimize for mean overall deltaG of dimers [True] or total dimer tally [False]? [Default: False] 

* **VERBOSE (-v –-verbose)** : Print all steps and swaps at the optimization step. [Default: False]

#### mw-assess-panel

* **PRIMERS (-i –-input)** : FASTA or CSV of primers. Sequence names must match the format .<#>.<FWD/REV> e.g., MACA01.0.FWD and MACA01.0.REV. If a CSV is provided, it must include 'PrimerID' and 'Sequence' fieldnames. 

* **ALL_DIMERS_dG (-a --alldimers_dg)** : Lower Gibbs free energy (deltaG) threshold for predicting non-end dimers. [Default: -8] 

* **END_DIMERS_dG (-e --enddimers_dg)** : deltaG threshold for predicting 3' end dimers. [Default: -4] 

* **BAD_DIMERS_dG (-b --baddimers_dg)** : deltaG threshold for counting dimers as particularly "bad". [Default: -10]

These are the basics but multiplex-wormhole has many more options- see [GitHub pages](https://mhallerud.github.io/multiplex_wormhole/multiplex-wormhole) for full functionality.


## Comments/Questions/Ideas
Please report any problems, questions, or potential enhancements in the [GitHub Issues](https://github.com/mhallerud/multiplex_wormhole/issues) page.
