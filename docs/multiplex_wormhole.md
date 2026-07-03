---
layout: default
nav_order: 1
parent: index
permalink: /multiplex-wormhole
---

# Multiplex Wormhole Pipeline
`multiplex-wormhole` is a wrapper function that performs the following submodules using defaults:

TEMPLATES (& KEEPLIST) --> [mw-primer-design](1_BatchPrimerDesign.md) --> [MFEprimer dimer](2_DimerPrediction.md) --> 
[mw-tabulate-dimers](3_TabulateDimers.md) --> [mw-mult-optimizations](5_MultipleRunOptimizations.md)

Note that `mw-mult-optimizations` also runs [mw-assess-panel](6_AssessPanel.md) at the end.


## Usage
### Command line syntax
```
multiplex-wormhole [-h] -t TEMPLATES -n NLOCI -o OUTDIR [-p PREFIX]
                   [-k KEEPLIST] [-r RUNS] [-i ITER] [-c CYCLES]
                   [-s SIMPLE] [-d] [-v]
```

### Python usage
```
# load module
import multiplex_wormhole as mw

# panel design example (showing defaults for optional parameters)
mw.multiplexWormhole(TEMPLATES="Input_Templates.csv", 
                     N_LOCI=50, 
                     OUTDIR="Test_MW", 
                     PREFIX="Test_MW_default",
                     KEEPLIST_FA="Keeplist.fa",
                     N_RUNS=10, ITERATIONS=1000, CYCLES=10, SIMPLE=5000, deltaG=False, VERBOSE=False)#optional
```

### Arguments
**TEMPLATES (-t --templates)** : Path to templates CSV. 
**NLOCI (-n --nloci)** : Final panel size (i.e., # primer pairs & # templates amplified).
**OUTDIR (-o --outdir)** : Filepath where output directory will be created and all outputs saved within a generated folder structure.
**PREFIX (-p --prefix)** : Prefix for all outputs. [Defaults to a timestamp if None provided]
**KEEPLIST_FA (-k --keeplist)** : Path to keeplist FASTA. [Default: None]
**N_RUNS (-r --runs)** : Number of optimization runs. [Default: 10]
**ITERATIONS (-i --iter)** : Number of simulated annealing iterations per cycle. [Default: 1000]
**CYCLES (-c --cycles) ** : Number of simulated annealing cycles per run. [Default: 10]
**SIMPLE (-s --simple)** : Number of simple iterative improvement iterations per run. [Default: 5000]
**deltaG (-d --deltaG)** : Optimize for mean overall deltaG of dimers [True] or total dimer tally [False]? [Default: False]
**VERBOSE (-v --verbose)** : Print all steps and swaps at the optimization step. [Default: False]


## Outputs
`multiplex-wormhole` sets up a directory structure in the designated `OUTDIR`, which contains the following outputs:

```
OUTDIR
├── 0_Inputs
│   ├── TEMPLATES
│   ├── KEEPLIST_FA
├── 1_PrimerDesign
│   ├── FilteredPrimers.csv
│   ├── FilteredPrimers.fa
│   ├── FilteredPrimers.log
│   ├── FilteredPrimers_plusKeeplist.fa
├── 2_PredictedDimers
│   ├── MFEprimerDimers.txt
│   ├── MFEprimerDimers_ends.txt
│   ├── PrimerPairInteractions_wide.csv
│   ├── PrimerPairInteractions_mean.csv
│   ├── PrimerPairInteractions_sum.csv
│   ├── PrimerPairInteractions_binary_sum.csv
│   ├── PrimerPairInteractions_binary_wide.csv
│   ├── PrimerPairInteractions.log
├── 3_OptimizedMultiplexes
    ├── PREFIX_RunSummary.csv
    ├── PREFIX_Runx_primers.fasta
    ├── Plots_Dimer_Load
    │   ├── PREFIX_Runx_DimerLoad.png
    ├── Final_Primers
    │   ├── PREFIX_Runx_primers.csv
    ├── Final_Dimers
    │   ├── PREFIX_Runx_dimers.csv
    ├── Trace_Dimer_Load
    │   ├── PREFIX_Runx_costsTrace.csv
    ├── Logfiles
        ├── PREFIX_Runx.log
```


[Next](1_BatchPrimerDesign.md)
