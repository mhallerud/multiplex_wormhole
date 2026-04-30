# multiplex_wormhole
*In silico* optimization for multiplex PCR assays that minimized predicted primer dimer loads. The pipeline was developed for genotyping by amplicon sequencing (i.e., reduced SNP panel) applications, however, the process is transferable to any multiplex PCR targeted sequencing approach. The impetus for multiplex_wormhole was genotyping noninvasive wildlife genetic samples. Default primer design settings are therefore conservative and tailored towards amplifying low concentration and degraded DNA such as that found in fecal and hair samples. 

Full documentation: [https://mhallerud.github.io/multiplex_wormhole](https://mhallerud.github.io/multiplex_wormhole/)

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
pixi add "python=3.9" #set python version
pixi shell #enter virtual env
```

### Install multiplex wormhole
```
pip install -i https://test.pypi.org/simple/ multiplex-wormhole
```
Note: Pixi/conda can be finicky... Dependening on your system, you may run into dependency errors here. If that happens, `exit` your virtual env and install the missing dependencies following the instructions below.

### Back-up installation option
If pip install doesn't work, you can also install manually by taking the following steps (from the command line):
1. Install Python dependencies to your virtual environment (replace "pixi add" with "conda install" if using conda):
```
pixi add primer3-py
pixi add pandas==1.4.4
pixi add numpy==1.24.4
pixi add matplotlib==3.5.2
```

Or, if not using a virtual environment:
```
pip install primer3-py==2.0.0
pip install pandas==1.4.4
pip install numpy==1.24.4
pip install matplotlib==3.5.2
```

2. Download source code from GitHub:
```
git clone https://github.com/mhallerud/multiplex_wormhole/
```

### Configuring the MFE primer binary
MFEprimer is used for dimer calculations. Multiplex wormhole is set up to automatically download and configure the binary file using the helpers/setup_mfeprimer.py script, take the following steps: Download the MFEprimer v3.2.7 version that fits your operating system [here](https://github.com/quwubin/MFEprimer-3.0/releases). Save the file to your multiplex_wormhole package directory (location can be found by running `pip show multiplex_wormhole`). Unzip the download (if zipped). Ensure the file can be executed by opening terminal or the command line in this directory and running `chmod +x mfeprimer*`.

Now you are ready to run multiplex wormhole!


## Full Workflow Example
See [multiplex_primer_design](multiplex_primer_design.py) to run the full pipeline step by step using the [example templates](examples/Input_Templates.csv).


## Comments/Questions/Ideas
Please report any problems, questions, or potential enhancements in the [GitHub Issues](https://github.com/mhallerud/multiplex_wormhole/issues) page.
