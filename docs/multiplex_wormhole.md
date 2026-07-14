---
layout: default
nav_order: 1
parent: index
permalink: /multiplex-wormhole
---

# Multiplex Wormhole Pipeline
`multiplex-wormhole` is a wrapper function that performs the following submodules using defaults:

[TEMPLATES](index.md#input-file-format) (& KEEPLIST) --> [mw-primer-design](1_BatchPrimerDesign.md) --> [MFEprimer dimer](2_DimerPrediction.md) --> 
[mw-tabulate-dimers](3_TabulateDimers.md) --> [mw-mult-optimizations](5_MultipleRunOptimization.md)

Note that `mw-mult-optimizations` also runs [mw-assess-panel](6_AssessPanel.md) at the end.


## Usage
### Command line syntax
```
multiplex-wormhole [-h] -t TEMPLATES -n NLOCI -o OUTDIR [-p PREFIX]
                   [-k KEEPLIST] [-r RUNS] [-i ITER] [-c CYCLES]
                   [-s SIMPLE] [--threads CPUs] [-d] [-v]
# arguments passed to primer design & dimer prediction steps:
                   [--tm-limit 45] [--dg-hairpins -2.0] [--dg-end-limit -4.0] [--dg-mid-limit -8.0]
                   [--primer3-settings '"{FWD_OVERHANG: "tcgtcggcagcgtcagatgtgtataagagacag"}'"] [--enable-broad] 
# arguments passed to optimization and plotting steps:
                   [--burnin 200] [--decay-rate 0.95] [--t-init None] [--t-final 0.01] [--prob_adj 2.0]
# arguments used in panel assessment:
                   [--dg-bad-limit -10]
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
                     # optional arguments passed to optimization:
                     N_RUNS=10, ITERATIONS=1000, CYCLES=10, SIMPLE=5000, deltaG=False, VERBOSE=False,
                     # primer design options                      
                     Tm_LIMIT=45, dG_HAIRPINS=-2, dG_END_LIMIT=-4,  dG_MID_LIMIT=-8, 
                     ENABLE_BROAD=False, PRIMER3_SETTINGS=None,
                     # multi-threading
                     THREADS=None,
                     # optimization / plotting options
                     BURNIN=200, DECAY_RATE=0.95, T_INIT=None, T_FINAL=0.01, PROB_ADJ=2,
                     # assessment limits
                     dG_BAD_LIMIT=-10)

```

### Primary Arguments
**TEMPLATES (-t --templates)** : Path to templates CSV. 
**NLOCI (-n --nloci)** : Final panel size (i.e., # primer pairs & # templates amplified).
**OUTDIR (-o --outdir)** : Filepath where output directory will be created and all outputs saved within a generated folder structure.
**PREFIX (-p --prefix)** : Prefix for all outputs. [Defaults to a timestamp if None provided]
**KEEPLIST_FA (-k --keeplist)** : Path to keeplist FASTA. [Default: None]
**N_RUNS (-r --runs)** : Number of optimization runs. [Default: 10]
**ITERATIONS (-i --iter)** : Number of simulated annealing iterations per cycle. [Default: 1000]
**CYCLES (-c --cycles) ** : Number of simulated annealing cycles per run. [Default: 10]
**SIMPLE (-s --simple)** : Number of simple iterative improvement iterations per run. [Default: 5000]
**THREADS (--threads)** : CPUs for multiprocessing- defaults to CPUs available on current machine. 
**deltaG (-d --deltaG)** : Optimize for mean overall deltaG of dimers [True] or total dimer tally [False]? [Default: False]
**VERBOSE (-v --verbose)** : Print all steps and swaps at the optimization step. [Default: False]

### Additional arguments passed to sub-modules:
**Tm_LIMIT (--tm-limit)** : Minimum melting temperature of predicted secondary structures, in Celsius, used in primer design. [Default: 45] 
**dG_HAIRPINS (--dg-hairpins)** : DeltaG threshold for discarding primers with hairpins, used in primer design. [Default: -2]
**dG_END_LIMIT (--dg-end-limit)** : DeltaG threshold for 3' end dimers used in primer design and dimer prediction. [Default: -4]
**dG_MID_LIMIT (--dg-mid-limit)** : DeltaG threshold for non-end dimers used in primer design and dimer prediction. [Default: -8]
**ENABLE_BROAD (--enable broad)** : Enable broader settings for primer design (if narrow settings fail for a given template)? [Default: False]
**PRIMER3_SETTINGS (--primer3-settings)** : Dictionary format of primer3 design settings. See [mw-primer-design](1_BatchPrimerDesign.md) for details.
**BURNIN (--burnin)** : Iterations used to sample cost space to set ASA temperatures for each cycle. [Default=200]
**DECAY_RATE (--decay-rate)** : Decay rate used to define ASA temperature schedule during optimization. [Default: 0.95]
**T_INIT (--t-init)** : Initial temperature used in adaptive simulated annealing algorithm. [Default: None (calculated adaptively)]
**T_FINAL (--t_final)** : Final temperature used in adaptive simulated annealing algorithm.[Default: 0.01]
**PROB_ADJ (--prob-adj)** : Multiplier used to adjust acceptance probabilities during adaptive simulated annealing. [Default: 2]
**dG_BAD_LIMIT (--dg-bad-limit)** : DeltaG threshold used to identify and count 'bad' dimers during panel assessment. [Default: -10]

## Outputs
`multiplex-wormhole` sets up a directory structure in the designated `OUTDIR`, which contains the following outputs:

```
OUTDIR
├── 0_Inputs : Copies of input files for this run
│   ├── TEMPLATES
│   ├── KEEPLIST_FA
├── 1_PrimerDesign : Outputs from mw-primer-design + mw-add-keeplist (if KEEPLIST provided)
│   ├── FilteredPrimers.csv : Designed primers with full details
│   ├── FilteredPrimers.fa : Sequences for filtered primers
│   ├── FilteredPrimers.log
│   ├── FilteredPrimers_plusKeeplist.fa
├── 2_PredictedDimers : Outputs from MFEprimer dimer and mw-tabulate-dimers
│   ├── MFEprimerDimers.txt : MFE primer dimer output
│   ├── MFEprimerDimers_ends.txt : MFE primer dimer output for 3' end dimers
│   ├── PrimerPairInteractions_wide.csv : pairwise dimer loads
│   ├── PrimerPairInteractions_mean.csv : mean deltaG primer pair across all others (if deltaG=True)
│   ├── PrimerPairInteractions_sum.csv : total dimer count per primer pair across all others (if deltaG=False)
│   ├── PrimerPairInteractions_binary_sum.csv : count of # interacting primer pairs per primer pair
│   ├── PrimerPairInteractions_binary_wide.csv : pairwise interactions between primer pairs (1=dimer, 0=no dimers)
│   ├── PrimerPairInteractions.log
├── 3_OptimizedMultiplexes : Outputs from mult-optimizations
    ├── PREFIX_RunSummary.csv : Dimer load summary for all optimization runs
    ├── PREFIX_Runx_primers.fasta : Primer sequences for "best" optimization runs
    ├── Plots_Dimer_Load : Dimer load trace plots across all optimization runs
    │   ├── PREFIX_Runx_DimerLoad.png
    ├── Final_Primers : CSVs of primer sequences for all optimization runs
    │   ├── PREFIX_Runx_primers.csv
    ├── Final_Dimers : Pairwise dimer tables for output panels from all optimization runs
    │   ├── PREFIX_Runx_dimers.csv
    ├── Trace_Dimer_Load : Dimer trace values for each iteration where a change was accepted during optimization runs
    │   ├── PREFIX_Runx_costsTrace.csv
    ├── Logfiles : Log files for all optimization runs
        └── PREFIX_Runx.log
```


[Next](1_BatchPrimerDesign.md)
