# Explore Optimization Parameters with `plot_SA_temps`

## Purpose
The optimization process relies on a simulated annealing algorithm defined by a variety of parameters. This script allows the user to explore the effects of these parameters on the algorithm's behavior when making risky decisions.

See Step 6_OptimizeMultiplexPrimerSet for more details on parameters and optimization with simulated annealing.

## Usage
### Python syntax
`import os`

`os.chdir('/multiplex_wormhole/`)

`from scripts.plot_SA_temps import main as plotSAtemps`

Checking default parameters for a fileset (provided values are defaults). These will be used to define MIN_DIMER and MAX_DIMER. Other parameters may be provided as desired:

`plotSAtemps(OUTPATH, PRIMER_FASTA, DIMER_SUMS=None, DIMER_TABLE, N_LOCI, WHITELIST=None, SEED=None, BURNIN=100)`

Testing parameters based on known cost values (provided values are defaults). If MIN_DIMER and MAX_DIMER are provided, filepaths from a known problem are ignored.

`plotSAtemps(OUTPATH, MIN_DIMER, MAX_DIMER, DECAY_RATE=0.95, T_INIT=None, T_FINAL=None, ADJUSTMENT=0.1)`

### Command line syntax
`cd /multiplex_wormhole/scripts`

`python3 plotSAtemps OUTPATH PRIMER_FASTA DIMER_SUMS DIMER_TABLE N_LOCI WHITELIST SEED BURNIN MIN_DIMER MAX_DIMER DECAY_RATE T_INIT T_FINAL ADJUSTMENT`

*Note: The command line requires all 13 arguments to be provided, however defaults such as 'None' can be provided where necessary. For example, to run a known set of parameters:*

`python3 plotSAtemps 'None' 'None' 'None' 'None' 'None' 'None' 'None' 'None' MIN_DIMER MAX_DIMER DECAY_RATE T_INIT T_FINAL ADJUSTMENT`

### Arguments
**OUTPATH** : Output directory path and prefix for output files.

**PRIMER_FASTA** : Filepath to FASTA file containing primer sequences to test. *Important: PrimerIDs in this file must match primer pair IDs in DIMER_SUMS and DIMER_TABLE!*

**DIMER_SUMS** : Filepath to CSV file containing dimer loads per primer pair. Output from tabulate_MFEprimer_dimer.

**DIMER_TABLE** : Filepath to CSV file containing pairwise primer dimer loads. Output from tabulate_MFEprimer_dimer.

**N_LOCI=None** : Number of loci amplified by multiplex primer set (including whitelist loci).

**WHITELIST** : Filepath to FASTA file containing primer sequences that MUST be included in final primer set. (Default: None)

**SEED** : Filepath to list of primers to use as the initial set of loci. This option overrides N_LOCI, so the number of loci in the SEED set will be the final number of loci in the primer set. (Default: None)

**MIN_DIMER** : Minimum 'bad' dimer change expected going from one iteration to the next.

**MAX_DIMER** : Maximum 'bad' dimer change expected moving from one iteration to the next.

**DECAY_RATE** : Variable between 0 and 1 to use in exponential decay function for temperature. Values close to 1 lead to slow temperature decay while decreasing values lead to more rapid temperature decay. (Default: 0.95)

**T_INIT=None** : Initial temperature where simulated annealing starts. If not provided, this value will be calculated from changes in dimer load. (Default: None)

**T_FINAL=None,** : Final temperature where simulated annealing ends. If not provided, this value will be calculated from changes in dimer load. (Default: None)

**BURNIN** : Number of iterations used to sample changes in dimer loads. Increase if dimer loads are highly variable, decrease if dimer loads are less variable. (Default: 100)

**ADJUSTMENT** : Value between 0 and 1 used as a multiplier in calculating T_INIT. Values closer to 1 result in high T_INIT and values closer to 0 result in T_INIT~T_FINAL. (Default: 0.1)


## Outputs
`OUTPATH`_TemperatureSchedule.png

![TestingASAparams_TemperatureSchedule](https://github.com/mhallerud/multiplex_wormhole/assets/43869036/1f31a6e3-67a8-48de-864c-a29477cfd5c7)

`OUTPATH`_DimerAcceptanceProbs.png

![TestingASAparams_DimerAcceptanceProbs](https://github.com/mhallerud/multiplex_wormhole/assets/43869036/c661c219-89ae-4724-89a1-60a63a059339)


[Previous](5_TabulateDimers.md)		[Next](6_OptimizeMultiplexPrimerSet.md)