# Tabulate Dimers with <tabulate_MFEprimer_dimers>

## Purpose
Converts text files output by MFEprimer dimer into tables of pairwise primer pair interactions and total secondary interactions per primer pair.
   
## Usage

## Outputs
- A N_PRIMERS x N_PRIMERS table which includes the total number of primer dimers estimated for all pairwise combinations of primer pairs.
- A N_PRIMERS x 2 table which includes the total number of primer dimers contributed by each primer pair.
- A N_PRIMERS X N_PRIMERS table which specifies pairwise primer pair interactions in binary (0 = no secondary structures between primer 1 and primer 2, 1 = at least 1 secondary structure between primer 1 and primer 2).
- A N_PRIMERS X 2 table which inclues the total number of primer pairs with which each primer pair interacts.

Optionally, tables can also be output for per-primer interactions between primers within the same pair if the <OUT_PRIMERPATH> argument is provided.

If you want to prioritize minimizing the number of pairwise interactions in a multiplex, use the binary outputs in the optimziation step. 
If you want to prioritize minimizing the total dimer load in a multiplex, use the count outputs in the optimization step.
