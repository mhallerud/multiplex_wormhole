---
title: Run Multiple Optimizations
layout: page
permalink: /run-multiple-optimizations
---
# Run Multiple Optimizations
The [optimization process](7_OptimizationProcess.md) includes randomness and outputs from runs will vary, even with equivalent input parameters. The variability in results will increase with increasing problem complexity, and therefore it's recommended to run the [optimization step](4_OptimizeMultiplexPrimerSet.md) multiple times and proceed with the best run. Variation between runs increases with problem complexity, so more runs are needed. This function automates running the optimization step multiple times for a given set of inputs.


## Usage
### Command line syntax
```
mult_optimizations -r RUNS -f PRIMER_FASTA -d DIMER_SUMS -t DIMER_TABLE -o OUTPATH -n NLOCI [-k KEEPLIST]
                    [--seed CSV] [-s GREEDY] [-i ITER] [-c CYCLES] [-b BURNIN]
                    [-decay-rate 2.0] [--t-init None] [--t-final 0.01] [--prob-adj 2.0]
                    [-g] [-v] [-m] [--threads 1]
                    [--dg-end-limit -4] [--dg-mid-limit -8] [--dg-bad-limit -10]
```

### Python syntax
```
import multiplex_wormhole as mw
# example with full inputs & their defaults
mw.multipleOptimizations(N_RUNS, PRIMER_FA, DIMER_SUMS, DIMER_TABLE, OUTPATH, N_LOCI, 
    deltaG=False, KEEPLIST=None, VERBOSE=False, SEED=None,
    GREEDY=5000, ITERATIONS=1000, CYCLES=10, BURNIN=100,
    DECAY_RATE=0.95, T_INIT=None, T_FINAL=None, PROB_ADJ=2)
```

### Parameters
* **PRIMER_FA (--primer-fasta -f)** : FASTA containing all primers that will be considered in the multiplex. *Important: Primer pair IDs in this file must match the IDs found in the DIMER_SUMS and DIMER_TABLE inputs!* These should also not include KEEPLIST loci, although the script will attempt to filter these out if found.
* **DIMER_SUMS (--dimer-sums -d)** : CSV containing dimer totals for each primer pair. 
* **DIMER_TABLE (--dimer-table -t)** : CSV containing pairwise dimer interactions between primer pairs. 
* **OUTPATH (--outpath -o)** : Prefix for output filepaths, including directory structure.
* **N_LOCI (--nloci -n)** : Size of the final multiplex panel, i.e., number of templates targeted.
* **KEEPLIST (--keeplist -k)** : FASTA containing primer pairs that MUST be included in the multiplex (e.g., primer pairs from a previous set, primer pairs for sex ID, etc.). *Important: These primer pairs must have been considered in dimer formation! Primer IDs must match IDs in the DIMER_SUMS and DIMER_TABLE inputs.* [Default: None]
* **deltaG (--deltaG -g)** : Minimize mean deltaG [True] or count of dimers- requires deltaG dimer tables! [Default: False]
* **SEED (--seed)** : CSV containing a multiplex primer set to use as a starting point in optimization, in the format output by the present optimization function. [Default: None]
* **GREEDY (--greedy -s)** : Number of iterations to run greedy local search algorithm. [Default: 5000]
* **ITERATIONS (--iter -i)** : Iterations per simulated annealing cycle, where all steps (accepted and rejected changes) are counted. [Default: 1000]
* **CYCLES (--cycles -c)** : Number of simulated annealing cycles per optimization run. [Default: 10]
* **BURNIN (--burnin -b)** : Number of samples taken of increased dimer costs used to calculate simulated annealing temperature schedule. Only steps that cause increased cost are counted so that this number equals the number of 'mistakes' sampled. [Default: 100]
* **DECAY_RATE (--decay-rate)** : Parameter for exponential decay function of simulated annealing temperatures. Values closer to 1 result in a slower decay from the initial to final temperature (and more "mistakes"), and values closer to 0.5 result in rapid decay towards the final temperature. [Default: 0.95]
* **T_INIT (--t-init)** : Initial temperature to use in fixed schedule simulated annealing. By default, T_INIT is adaptively set based on the problem at hand. Higher initial temperatures means that more of the cost optimization space is explored, but more "mistakes" will also be allowed in the process. (default: None) By default, T_FINAL is set adaptively based on the problem at hand where: `T_FINAL=MIN_DIMERS + DIMER_ADJ * (MAX_DIMERS - MIN_DIMERS)` with MAX_DIMERS and MIN_DIMERS calculated from changes observed during the BURNIN stage. [Default: None -calculated from data]
* **T_FINAL (--t-final)** : Final temperature to use in fixed schedule simulated annealing. As temperatures approach 0, simulated annealing allows fewer 'mistakes' and converges with a greedy algorithm. [Default: 0.1]
* **PROB_ADJ (--prob-adj)** : Multiplier used to adjust dimer acceptance probabilities. Increased values result in lower dimer acceptance probabilities at the cost of exploring less of the cost optimization space. [Default: 2]
* **VERBOSE (-v)** : Print updates as algorithm proceeds? [Default: False]
* **MAKEPLOT (-m)** : Make simulated annealing temperature schedule plots? Runs [explore optimization parameters](4A_ExploreOptimParameters.md) function. [Default: False]
* **THREADS (--threads)** : Number of processors to use for multiprocessing. [Default: None]
* **dG_END_LIMIT (--dg-end-limit)** : DeltaG threshold used to count 3' end dimers during panel assessment. [Default: -4]
* **dG_MID_LIMIT (--dg-mid-limit)** : DeltaG threshold used to count non-end dimers during panel assessment. [Default: -8]
* **dG_BAD_LIMIT (--dg-bad-limit)** : DeltaG threshold used to count 'bad' dimers during panel assessment. [Default: -10]


## Outputs
Outputs are copied into directories with prefixes as `OUTPATH`_Run`N`...:
* **Final_Primers** : CSV of primer sequences in optimized multiplex.
* **Final_Dimers** : CSV of pairwise dimer loads in the optimized multiplex.
* **Trace_Dimer_Load** : CSV of dimer load reported at each swap in optimization
* **Dimer_Load_Plots** : Dimer load plots
* **Logfiles** : stdout log for each run.

The dimer load per optimized multiplex is output to `OUTPATH`_RunSummary.csv. Then, the best set(s) are output to FASTAs in the output directory as Run`N`_`OUTPATH`_primers.fa


[Previous](4_OptimizeMultiplexPrimerSet.md)           [Next](6_AssessPanel.md)
