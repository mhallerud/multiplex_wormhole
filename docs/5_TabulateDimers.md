# Tabulate Dimers with `tabulate_MFEprimer_dimers`

## Purpose
Converts text files output by MFEprimer dimer into tables of pairwise primer pair interactions and total secondary interactions per primer pair.
   
## Usage
### Python syntax
`import os`

`os.chdir('/multiplex_wormhole')`

`from scripts.tabulate_MFEprimer_dimers import main as tabulateDimers`

`tabulateDimers(ALL_DIMERS, END_DIMERS, OUTPATH, OUTPRIMERPATH="False")`

### Command line syntax
`cd multiplex_wormhole/scripts`

`python3 tabulate_MFEprimer_dimers.py ALL_DIMERS END_DIMERS OUTPATH False`

### Arguments
**ALL_DIMERS** : Text file output for all dimers output from MFE primer.

**END_DIMERS** : Text file output for end dimers output from MFE primer. NOTE: If you do not have any end dimers, you can provide a blank text file here.

**OUTPATH** : Directory path and prefix for output tables of primer pair interactions.

**OUTPRIMERPATH** (optional): Directory path and prefix for output tables of primer interactions.


## Outputs
1. A N_PRIMERS x N_PRIMERS table which includes the total number of primer dimers estimated for all pairwise combinations of primer pairs.
2. A N_PRIMERS x 2 table which includes the total number of primer dimers contributed by each primer pair.
3. A N_PRIMERS X N_PRIMERS table which specifies pairwise primer pair interactions in binary (0 = no secondary structures between primer 1 and primer 2, 1 = at least 1 secondary structure between primer 1 and primer 2).
4. A N_PRIMERS X 2 table which inclues the total number of primer pairs with which each primer pair interacts.

Optionally, tables can also be output for per-primer interactions between individual primers (not just pairs) if the `OUT_PRIMERPATH` argument is provided.

If you want to prioritize minimizing the number of pairwise interactions in a multiplex, use the binary outputs in the optimziation step. 
If you want to prioritize minimizing the total dimer load in a multiplex, use the count outputs in the optimization step.


[Previous](4_DimerPrediction.md)		[Next](6A_ExploreOptimParameters.md)