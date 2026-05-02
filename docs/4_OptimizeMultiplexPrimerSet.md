---
title: "Optimize Multiplex Primer Set"
layout: default
permalink: /optimize-multiplex-pcr
nav_order: 1
parent: index
---
# Optimize Multiplex PCR Primers
A panel that amplifies "N" target loci is selected by minimizing off-target interactions among primer pairs. Optionally, a "keeplist" set of primers (containing primers from a previous assay, primers for important loci, etc.) can be provided that *must* be included in the final primer set. 

[See details on the optimization process](7_OptimizationProcess.md) to better understand parameter settings.


## Usage
### Python syntax
```
import multiplex_wormhole as mw

# with minimal inputs:
mw.optimizeMultiplex(PRIMER_FASTA, DIMER_SUMS, DIMER_TABLE, OUTPATH, N_LOCI, deltaG=False)

# with all parameters & their defaults:
mw.optimizeMultiplex(PRIMER_FASTA, DIMER_SUMS, DIMER_TABLE, OUTPATH, N_LOCI, KEEPLIST=None, deltaG=False, SEED=None, #data args
SIMPLE=5000, ITERATIONS=1000, CYCLES=10, BURNIN=100, DECAY_RATE=0.95, T_INIT=None, T_FINAL=None, PROB_ADJ=2, MAKEPLOT=False,  VERBOSE=False, RNG=12345) #sim anneal params
```

### Command line syntax
```
cd ~/multiplex_wormhole #navigate to scripts
python3 optimize_multiplex.py -f PRIMER_FASTA -d DIMER_SUMS -t DIMER_TABLE -o OUTPATH -n NLOCI [-k KEEPLIST] [-e SEED] [-s SIMPLE] /
[-i ITER] [-c CYCLES] [-b BURNIN] [-r DECAY_RATE] [-x TEMP_INIT] [-l TEMP_FINAL] [-a PROB_ADJ] [-g] [-v] [-m]
```

### Parameters
**PRIMER_FASTA (-f)** : FASTA containing all primers that will be considered in the multiplex. *Important: Primer pair IDs in this file must match the IDs found in the DIMER_SUMS and DIMER_TABLE inputs!* These should also not include KEEPLIST loci, although the script will attempt to filter these out if found.

**DIMER_SUMS (-d)** : CSV containing dimer totals for each primer pair. 

**DIMER_TABLE (-t)** : CSV containing pairwise dimer interactions between primer pairs. 

**OUTPATH (-o)** : Prefix for output filepaths, including directory structure.

**N_LOCI (-n)** : Size of the final multiplex panel, i.e., number of templates targeted.

**KEEPLIST (-k)** : FASTA containing primer pairs that MUST be included in the multiplex (e.g., primer pairs from a previous set, primer pairs for sex ID, etc.). *Important: These primer pairs must have been considered in dimer formation! Primer IDs must match IDs in the DIMER_SUMS and DIMER_TABLE inputs.* [Default: None]

**deltaG (-g)** : Minimize mean deltaG [True] or count of dimers- requires deltaG dimer tables! [Default: False]

**SEED (-e)** : CSV containing a multiplex primer set to use as a starting point in optimization, in the format output by the present optimization function. [Default: None]

**SIMPLE (-s)** : Number of iterations to run simple iterative improvement optimization. [Default: 5000]

**ITERATIONS (-i)** : Number of iterations to run per simulated annealing cycle, where all steps (accepted and rejected changes) are counted. [Default: 1000]

**CYCLES (-c)** : Number of simulated annealing cycles to run. [Default: 10]

**BURNIN (-b)** : Number of samples taken of increased dimer costs used to calculate simulated annealing temperature schedule. Only steps that cause increased cost are counted so that this number equals the number of 'mistakes' sampled. [Default: 100]

**DECAY_RATE (-r)** : Parameter for exponential decay function of simulated annealing temperatures. Values closer to 1 result in a slower decay from the initial to final temperature (and more "mistakes"), and values closer to 0.5 result in rapid decay towards the final temperature. [default: 0.95]

**T_INIT (-x)** : Initial temperature to use in fixed schedule simulated annealing. By default, T_INIT is adaptively set based on the problem at hand. Higher initial temperatures means that more of the cost optimization space is explored, but more "mistakes" will also be allowed in the process. (default: None) By default, T_FINAL is set adaptively based on the problem at hand where: `T_FINAL=MIN_DIMERS + DIMER_ADJ * (MAX_DIMERS - MIN_DIMERS)` with MAX_DIMERS and MIN_DIMERS calculated from changes observed during the BURNIN stage. [Default: None -calculated from data]

**T_FINAL (-l)** : Final temperature to use in fixed schedule simulated annealing. As temperatures approach 0, simulated annealing allows fewer 'mistakes' and converges with simple iterative improvement. [Default: 0.1]

**PROB_ADJ (-a)** : Multiplier used to adjust dimer acceptance probabilities. Increased values result in lower dimer acceptance probabilities at the cost of exploring less of the cost optimization space. [Default: 2]

**VERBOSE (-v)** : Print updates as algorithm proceeds? [Default: False]

**MAKEPLOT (-m)** : Make simulated annealing temperature schedule plots? Runs [explore optimization parameters](4A_ExploreOptimParameters.md) function. [Default: False]

**RNG*** : Random number generator used to ensure reproducibility.


## Outputs
1. **`OUTPATH`_dimers.csv** : Table with pairwise dimer loads of primer pairs in optimized multiplex.
2. **`OUTPATH`_primers.csv** : Table with Primer IDs and (adapter-ligated) sequences. 
3. **`OUTPATH`_costsTrace.csv** : Trace of dimer load cost, simulated annealing temperature, and iterations as the algorithm progressed. Only accepted changes are reported.
4. **`OUTPATH`_DimerLoad.png** : Plot of dimer load trace across iterations. See example and details [here](https://mhallerud.github.io/multiplex_wormhole/optimization-process#understanding-the-dimer-trace-plot).



[Previous](4A_ExploreOptimParameters.md)
