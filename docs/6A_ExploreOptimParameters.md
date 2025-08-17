# Explore Optimization Parameters with `plot_SA_temps`

## Purpose
The optimization process relies on a simulated annealing algorithm defined by a variety of parameters. This script allows the user to explore the effects of these parameters on the algorithm's behavior when making risky decisions.

See [Step 6_OptimizeMultiplexPrimerSet](6_OptimizeMultiplexPrimerSet.md) for more details on parameters and optimization with simulated annealing. This function is also called within the optimization function to plot parameters used in simulated annealing.

## Usage
### Python syntax
```
import os
os.chdir('/multiplex_wormhole/`)
from scripts.plot_SA_temps import main as plotSAtemps

# Checking default parameters for a fileset (provided values are defaults). These will be used to define MIN_DIMER and MAX_DIMER. Other parameters may be provided as desired:
plotSAtemps(OUTPATH, PRIMER_FASTA, DIMER_SUMS=None, DIMER_TABLE, N_LOCI, KEEPLIST=None, SEED=None, BURNIN=100)

# Testing parameters based on known cost values (provided values are defaults). If MIN_DIMER and MAX_DIMER are provided, filepaths from a known problem are ignored.
plotSAtemps(OUTPATH, MIN_DIMER=1, MAX_DIMER=5, DECAY_RATE=0.95, T_INIT=None, T_FINAL=None, DIMER_ADJ=0.1, PROB_ADJ=1.5)
```

### Command line syntax
```
cd /multiplex_wormhole/scripts
python3 plotSAtemps OUTPATH PRIMER_FASTA DIMER_SUMS DIMER_TABLE N_LOCI KEEPLIST SEED MIN_DIMER MAX_DIMER DECAY_RATE T_INIT T_FINAL BURNIN DIMER_ADJ PROB_ADJ

#Note: The command line requires all 13 arguments to be provided, however defaults such as 'None' can be provided where necessary. For example, to run a known set of parameters:*
python3 plotSAtemps OUTPATH 'None' 'None' 'None' 'None' 'None' 'None' MIN_DIMER MAX_DIMER DECAY_RATE T_INIT T_FINAL BURNIN DIMER_ADJ PROB_ADJ
```

### Arguments
**OUTPATH** : Output directory path and prefix for output files.

**PRIMER_FASTA** : Filepath to FASTA file containing primer sequences to test. *Important: PrimerIDs in this file must match primer pair IDs in DIMER_SUMS and DIMER_TABLE!* (Default: None)

**DIMER_SUMS** : Filepath to CSV file containing dimer loads per primer pair. Output from tabulate_MFEprimer_dimer. (Default: None)

**DIMER_TABLE** : Filepath to CSV file containing pairwise primer dimer loads. Output from tabulate_MFEprimer_dimer. (Default: None)

**N_LOCI** : Target panel size (includes keeplist loci), i.e. the number of unique amplicons resulting from the multiplex primer set. (Default: None)

**KEEPLIST** : Filepath to FASTA file containing primer sequences that MUST be included in final primer set. (Default: None)

**SEED** : CSV to primer set used as initial loci, in format of multiplex_wormhole output *_primers.csv. This option overrides N_LOCI, so the number of loci in the SEED set will be the final number of loci in the primer set. (Default: None)

**MIN_DIMER** : Minimum 'bad' dimer change expected going from one iteration to the next. (Default: None)
Generally MIN_DIMER=1

**MAX_DIMER** : Maximum 'bad' dimer change expected moving from one iteration to the next. (Default: None)

**DECAY_RATE** : Variable between 0 and 1 to use in exponential decay function for temperature. Values close to 1 lead to slow temperature decay while decreasing values lead to more rapid temperature decay. (Default: 0.98)

**T_INIT** : Initial temperature where simulated annealing starts. If not provided, this value will be calculated from changes in dimer load. (Default: None- set adaptively)
Adaptively calculated as `MIN_DIMER + DIMER_ADJ * (MAX_DIMER - MIN_DIMER)`

**T_FINAL** : Final temperature where simulated annealing ends. If not provided, this value will be calculated from changes in dimer load. (Default: 0.1)

**BURNIN** : Number of iterations used to sample changes in dimer loads. Increase if dimer loads are highly variable, decrease if dimer loads are less variable. (Default: 100)

**DIMER_ADJ** : Value between 0 and 1 used as a multiplier in calculating T_INIT. Values closer to 1 result in high T_INIT and values closer to 0 result in T_INIT~T_FINAL. (Default: 0.1)

**PROB_ADJ** : Multiplier to adjust dimer acceptance probability, where 1=no adjustment and higher values reduce dimer acceptance probability. (Default: 2) 
The probability of accepting an increase in dimer load follows `e^(-PROB_ADJ*cost / temperature)`

## Outputs
`OUTPATH`_TemperatureSchedule.png

![TestingASAparams_TemperatureSchedule](https://github.com/mhallerud/multiplex_wormhole/assets/43869036/1f31a6e3-67a8-48de-864c-a29477cfd5c7)

`OUTPATH`_DimerAcceptanceByTemp.png

![TestingASAparams_DimerAcceptanceProbs](https://github.com/mhallerud/multiplex_wormhole/assets/43869036/c661c219-89ae-4724-89a1-60a63a059339)

`OUTPATH`_DimerAcceptanceByIteration.png


[Previous](5_TabulateDimers.md)		[Next](6_OptimizeMultiplexPrimerSet.md)