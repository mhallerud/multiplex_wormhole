---
title: Multiplex Optimization Process
layout: page
permalink: optimization-process
---

# The Optimization Process

[basic idea]


## Cost Functions
multiplex_wormhole can describe "dimer load" using two distinct cost functions: 

With the **standard approach** (`deltaG=False`), cost of the multiplex is calculated as the total count of  dimers between primer pairs within the set. The cost table for this approach can be either the count of dimers per pairwise interaction (i.e., to minimize for *total* dimer load) or a binary pairwise table indicating which primer pairs interact (i.e., to minimize the *number of interacting pairs*). By default with this approach, multiplex wormhole minimizes total dimer load.

Alternatively, with the **deltaG approach** (`deltaG=True`), cost of the multiplex reflects the "badness" or strength of dimers within the multiplex. For each combination of primer pairs, the minimum deltaG (i.e., worst or strongest) is identified. Then, this is averaged across pairwise combinations within the multiplex. The maximum "mean deltaG" (i.e., closest to zero - reflecting weakest dimers) is then optimized for. 

For simple problems, both approaches will return a final predicted dimer load of zero. However, for complex problems where dimerization is inevitable, the deltaG option allows strength of dimers to be considered rather than treating all dimers equally. This may be particularly valuable when building large multiplexes. 


## Optimization Procedure
## 1. Selecting an initial primer set: Pseudo-greedy algorithm
1. The best primer pairs (i.e., the primer pairs with the lowest dimer load across all other candidate primer pairs) are identified for each template.

2. An initial set of `N_LOCI` primer pairs is then selected by choosing the best (i.e., lowest cumulative dimer load) templates from the template-specific "bests". If a keeplist is provided, the keeplist primer pairs are included in the initial set and `N_LOCI`-`N_KEEPLIST_PAIRS` is used to fill out the set.


## 2. Exploring the optimization space: Simulated annealing algorithm
*Simulated annealing* is an optimization procedure that allows the algorithm to make "mistakes" as it proceeds. Short-term "mistakes" (i.e., changes that result in a worst set) may make the primer set worst at first, but help to overcome local minima in the optimization space using "hill-climbing" behavior. Simulated annealing proceeds along a schedule where riskier decisions (i.e., decisions that result in worse outcomes from the previous iterations) are allowed in early iterations, but as the algorithm proceeds, riskier decisions are allowed with lower probability. 

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

## 3. Fine-tuning the optimized set: Simple iterative improvement algorithm
The best set (lowest dimer load) from simulated annealing is used as input into an iterative improvement algorithm. In iterative improvement, the worst primer pair is identified and swapped with a random alternative. If the new primer pair has fewer dimers than the previous, then it is kept and the new primer set is used as the starting point for the next iteration. If not, all other alternatives for this worst primer pair are tried until either a better option is found or alternatives are exhausted. Then, the algorithm proceeds to the next worst primer pair. This process continues until either A) the primer set is fully optimized (0 dimers), B) no alternatives can improve the set anymore, or C) the user-defined number of iterations `SIMPLE` have passed.


## Additional Tips & Tricks
* See the `plot_SA_temps` script to test the effects of `DECAY_RATE`, `T_INIT`, `T_FIXED`, `DIMER_ADJ`, and `PROB_ADJ` on the probability of accepting a range of dimer values as the simulated annealing algorithm proceeds.
* The optimization process is inherently random and results will vary across runs. As such, it's recommended to run multiple times and select the "best" result. `multiple_run_optimization` provides an easy wrapper for running and summarizing across multiple runs with standard inputs. Variability in results will increase with problem complexity, and more runs are recommended for very complex problems.
* For simple problems (e.g., selecting 50 primer pairs from 200 candidates), the iterative improvement is sufficient and more efficient. To run the simple iterative improvement algorithm, set `SIMPLE` high and `ITERATIONS=0` to bypass simulated annealing. 
* For challenging problems (e.g., many options and/or many `N_LOCI`), you may consider successive optimization runs. The `SEED` argument in the `optimizeMultiplex` function is provided to allow an output from a previous run to be used as the "initial panel" for a new optimization run. This is not the same as just increasing the number of iterations because the simulated annealing temperature schedule will be recalculated for this set and simulated annealing will possibly allow local minima to be overcome.
* When panels are designed to provide multiple levels of information (e.g., individual ID + sex + hybrid status), multiple optimization runs will be necessary to get the desired information content for each inferential objective. 
* Optimization is constrained by the number of options available relative to the number of primer pairs desired. The number of candidates required to design a panel with minimal dimer load drastically increases as desired panel size increases.
