---
title: Multiplex Optimization Process
layout: page
permalink: optimization-process
---

# The Optimization Process

An initial primer set is selected using a pseudo-greedy algorithm where the primer pairs with the cumulative lowest dimer load (across all loci provided) are selected, then adaptive simulated annealing is used to explore the optimization space around this initial primer set by randomly swapping out primer pairs and keeping improvements while allowing for 'mistakes' that may improve the primer set in the long run, and finally the best primer set found during adaptive simulated annealing is entered into a simple iterative improvement algorithm where the worst loci are swapped for better alternatives.

### 1. Selecting an initial primer set: Pseudo-greedy algorithm
A. The best primer pairs (i.e., the primer pairs with the lowest dimer load) are identified for each locus.
*NOTE: Dimer load can be calculated either as the total number of primer dimers that loci contributes, or as the number of binary primer pair interactions (1 = primer pair has at least one primer dimer, 0 = primer pair has no interactions) the locus contributes. By default, multiplex_wormhole considers the binary primer pair interactions. To optimize based on total primer dimers instead, use the *_PrimerPair_interactions_sum.csv and *_PrimerPair_interactions_wide.csv as the DIMER_SUMS and DIMER_TABLE inputs.* 

B. An initial set of primer pairs is selected by choosing the best (i.e., lowest dimer load) `N_LOCI` primer pairs from the locus-specific bests.

### 2. Exploring the optimization space: Adaptive simulated annealing algorithm
*Simulated annealing* is an optimization procedure that allows the algorithm to make "mistakes" as it proceeds. Short-term "mistakes" (i.e., changes that result in a worst set) may make the primer set worst at first but help to overcome local minimums in the optimization space. Simulated annealing proceeds along a schedule where riskier decisions (i.e., decisions that result in worse outcomes) are allowed in early iterations, but as the algorithm proceeds riskier decisions are allowed with lower probability since local optima should have been overcome at this point. 

At each iteration in the algorithm, a primer pair is randomly swapped with an alternative. The total dimer load (i.e., cost) of the new primer set is calculated and compared to the cost of the previous set. Based on the simulated annealing temperature and the calculated cost, the change will be accepted and the algorithm will use this primer set moving forward. All changes that improve the primer set (i.e., lower dimer cost) are accepted, while the probability of accepting a change that results in a worse set is governed by a temperature schedule. This process is repeated 1000s of times (based on `ITERATIONS`), and as the algorithm proceeds the primer set should steadily improve. More iterations will be needed for more challenging problems.

In multiplex_wormhole, an adaptive simulated annealing approach is used to set the temperature schedule. After sampling `BURNIN` "mistakes" where a random change resulted in a worse primer set, temperatures are set based on the least bad (`MIN_DIMERS`) and worst (`MAX_DIMERS`) changes observed. Specifically,

`T_INIT = MIN_DIMERS + DIMER_ADJ * (MAX_DIMERS - MIN_DIMERS)` with 0.1 as the default `DIMER_ADJ` parameter.

`T_FINAL = 0.1`

Higher temperatures correspond with allowing riskier decisions, and lower temperatures correspond with less risky decisions. A negative exponential function is used for the temperature gradient so that temperatures decline rapidly and most of the algorithm's time is spent in the "less risky" optimization space. The temperature schedule follows the function:

`Temperature = (T_INIT-T_FINAL) * 0.95^step + T_FINAL` with 0.95 as the default `DECAY_RATE` parameter.

Temperatures proceed to the next 'step' based on the `PARTITIONS` parameter which defines how many temperature changes will occur during simulated annealing. By default, with 1000 `PARTITIONS` and 5000 `ITERATIONS`, the temperature is reduced every 5 iterations.

At each iteration, the probability of accepting a change that increases dimer load is governed by the formula:

`e^(-PROB_ADJ*cost / temperature)`

where PROB_ADJ is set by the user, temperature is the simulated annealing temperature at the current step, and cost is the increase in dimer load relative to the previous iteration. 

### 3. Fine-tuning the optimized set: Simple iterative improvement algorithm
The best set (lowest dimer load) from simulated annealing is used as input into an iterative improvement algorithm. In iterative improvement, the worst primer pair is identified and swapped with a random alternative. If the new primer pair has fewer dimers than the previous, then it is kept and the new primer set is used as the starting point for the next iteration. If not, all other alternatives for this worst primer pair are tried until either a better option is found or alternatives are exhausted. Then, the algorithm proceeds to the next worst primer pair. This process continues until either A) the primer set is fully optimized (0 dimers), B) no alternatives can improve the set anymore, or C) the user-defined number of iterations `SIMPLE` have passed.


## Additional notes and recommendations on optimization...
See the `plot_SA_temps` script to test the effects of `DECAY_RATE`, `T_INIT`, `T_FIXED`, `DIMER_ADJ`, and `PROB_ADJ` on the probability of accepting a range of dimer values as the simulated annealing algorithm proceeds.

The optimization process is inherently random and results will vary across runs. As such, it's recommended to run multiple times and select the "best" result. `multiple_run_optimization` provides an easy wrapper for running `optimize_primers` multiple times in Python, however if you have cluster access I recommend running many runs on separate nodes instead.

For simple problems (e.g., selecting 50 primer pairs from a 200-locus set), the iterative improvement is sufficient and more efficient. To run the simple iterative improvement algorithm, set `SIMPLE` high and `ITERATIONS=0` to bypass simulated annealing. 

For challenging problems (e.g., many options and/or many `N_LOCI`), you may consider proceeding stepwise where first, an initial small set of loci is optimized from a subset of available loci, then next provided as a `KEEPLIST` into a run where `N_LOCI` is increased and the full set of options is provided.

Optimization is constrained by the number of options available relative to the number of primer pairs desired. In general, a ratio of 1:4 for desired:wanted primer pairs seems to work best. Providing more options will slow the process down, while providing too few options may result in unwanted high dimer loads. Providing an overabundance of options may also saturate the optimization space with an exponential number of extra options, thereby making the algorithm not only inefficient but also ineffective and preventing a good set from being found. For complex problems where dimer loads are consistent across all loci, more options are necessary. For problems where pairwise dimer loads vary widely, fewer options may be required. 

There is a tradeoff between computation time in simulated annealing (`ITERATIONS`) and optimization. For complex problems, longer runtimes and more runs will help explore more of the space. When in doubt, try plotting the trace of costs. If it is continuing to decrease at the tail end of iterations, increase the iterations. If the final cost is very variable between runs, increase the number of runs. You can also try subsetting the number of options to simplify the problem. Finally, if none of those work, you may have to consider aiming for a smaller target set.

