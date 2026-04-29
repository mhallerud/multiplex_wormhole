# multiplex_wormhole
*In silico* optimization for multiplex PCR assays that minimized predicted primer dimer loads. The pipeline was developed for genotyping by amplicon sequencing (i.e., reduced SNP panel) applications, however, the process is transferable to any multiplex PCR targeted sequencing approach. The impetus for multiplex_wormhole was genotyping noninvasive wildlife genetic samples. Default primer design settings are therefore conservative and tailored towards amplifying low concentration and degraded DNA such as that found in fecal and hair samples. 

Documentation webpage: [https://mhallerud.github.io/multiplex_wormhole](https://mhallerud.github.io/multiplex_wormhole/)

## Installation
### Set up a virtual environment
multiplex_wormhole was built and tested on MacOS with Python v3.9.13 in the Spyder IDE managed under Anaconda-Navigator. For those new to Python or with existing python packages, [Anaconda](https://www.anaconda.com/products/navigator) is the recommended virtual environment manager. For a conda virtual environment within your working directory:
```
conda create -n py39 python=3.9 #create new virtual env w/ python v3.9
conda activate py39 #enter virtual env
```
Some clusters used pixi instead of conda environments:
```
pixi init #initialize virtual env
pixi add "python==3.9" #set python version 
pixi shell #enter virtual env
```

### Setting up dependencies
The following Python packages are required and can be installed after entering your virtual env:
```
pip install primer3-py==2.0.0
pip install pandas==1.4.4
pip install numpy==1.24.4
pip install matplotlib==3.5.2
```
General Python modules required (these come pre-installed with most Python installations): os, sys, csv, string, subprocess, random, math, signal, gc, itertools, shutil, glob, datetime, argparse

### Download multiplex wormhole
From the command line:
```
git clone https://github.com/mhallerud/multiplex_wormhole/
```

### Configure up MFE primer binaries
MFEprimer is used for dimer calculations and can be configured with the setup_mfeprimer.py script. In python:
```
# load the function
import sys
sys.path.append("pathtoyour/multiplex_wormhole/src")
from setup_mfeprimer import main as setup_mfeprimer
# run
setup_mfeprimer()
```
From the command line:
```
# navigate into your mw directory
cd pathtoyour/multiplex_wormhole/src
python3 setup_mfeprimer.py 
```

If this fails, you can also manually download [here](https://github.com/quwubin/MFEprimer-3.0/releases). Then unzip and move the binaries into your /multiplex_wormhole/src directory. MFEprimer Version 3.2.7 was used during development. 


## Full Workflow Example
See [multiplex_primer_design](multiplex_primer_design.py) to run the full pipeline step by step using the [example templates](examples/Input_Templates.csv).


## Comments/Questions/Ideas
Please report any problems or potential enhancements in the GitHub Issues page.
