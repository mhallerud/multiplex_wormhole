# Optimize a Primer Set for Multiplex PCR with <optimize_primers>

## Purpose
A set of primers for "N" loci is selected that minimizes the number of negative interactions between primer pairs. Optionally, a "whitelist" set of primers (containing primers from a previous assay, primers for important loci, etc.) may be provided that *must* be included in the final primer set.

## The optimization process
An initial primer set is selected using a pseudo-greedy algorithm where the primer pairs with the cumulative lowest dimer load (across all loci provided) are selected, then adaptive simulated annealing is used to explore the optimization space around this initial primer set by randomly swapping out primer pairs and keeping improvements while allowing for 'mistakes' that may improve the primer set in the long run, and finally the best primer set found during adaptive simulated annealing is entered into a simple iterative improvement algorithm where the worst loci are swapped for better alternatives.

### Selecting an initial primer set: Pseudo-greedy algorithm

### Exploring the optimization space: Adaptive simulated annealing algorithm

### Fine-tuning optimized set: Simple iterative improvement algorithm
- The best primer pairs (i.e., the primer pairs with the lowest dimer load) are identified for each locus.
  NOTE: Dimer load can be calculated either as the total number of primer dimers that loci contributes, or as the number of binary primer pair interactions (1 = primer pair has at least one primer dimer, 0 = primer pair has no interactions) the loci contributes. By default, the algorithm considers the binary primer pair interactions. To optimize based on total primer dimers instead, use the *_PrimerPair_interactions_sum.csv and *_PrimerPair_interactions_wide.csv as the DIMER_SUMS and DIMER_TABLE inputs. 
- An initial set of primer pairs is selected by choosing the best (i.e., lowest dimer load) "N" primer pairs from the locus-specific bests.
- The algorithm then iteratively improves this initial set of primer pairs by switching the worst (i.e., highest dimer load) primer pair with a random alternative (either for the same locus, or from a locus that is not currently included in the set).
- If the change resulted in an improvement (i.e., fewer TOTAL dimer load for the set), then this new set is used as the starting point for the next iteration.
- If the change did not improve the set, then other alternatives for this worst pair are tested until a replacement that improves the set is found. If no alternatives improve the set, then this primer pair is blacklisted (i.e., kept in the set but a message will be returned suggesting to possibly remove this pair as an option) and the algorithm moves on to the next worst primer pair.
- The algorithm stops when one of the following conditions are met:
- The total primer dimer load is 0, and therefore there is no longer a need to optimize.
- No further improvements can be made to the set (i.e., all remaining primer pairs have been blacklisted- no alternatives to these will result in improvements).
- A user-defined number of iterations is reached.

### Default behavior
BURNIN = 100
ADJUSTMENT = 0.1
PARTITIONS = 1000
ITERATIONS = 5000
SIMPLE = 500
DECAY_RATE = 0.95

Temperatures:
T_init = min(change) + ADJUSTMENT * (max(change) - min(change))
T_end = min(change)
Temp = (T_init-T_end)*DECAY_RATE**step+T_end
    if ITERATIONS > PARTITIONS:
        Tspace = (T_init - T_end)/PARTITIONS
        T_iter = Tspace*ITERATIONS
    else:
        Tspace = (T_init - T_end)/ITERATIONS
        T_iter = 1    

