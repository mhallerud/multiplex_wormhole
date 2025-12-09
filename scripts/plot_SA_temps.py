#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: Plot Simulated Annealing Temperatures
Purpose: Assessing effects of simulated annealing parameter choices on acceptance
probability of dimer loads.

Created on Sun Mar 17 13:21:40 2024
@author: maggiehallerud
"""


# load dependencies
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import random as rand
import csv
import sys
import math

sys.path.append(os.path.dirname(__file__))
from optimize_primers import *



def main(OUTPATH, PRIMER_FASTA=None, DIMER_SUMS=None, DIMER_TABLE=None, N_LOCI=None, KEEPLIST=None, SEED=None, 
         MIN_DIMER=None, MAX_DIMER=None, DECAY_RATE=0.98, T_INIT=None, T_FINAL=0.1, BURNIN=100, DIMER_ADJ=0.1, PROB_ADJ=2):
    """
    PRIMER_FASTA : Fasta path
        Contains primer IDs and sequences
    DIMER_SUMS : CSV path
        Sum of interactions per primer pair (binary interactions recommended)
    DIMER_TABLE : CSV path
        Pairwise interaction table of primer pairs (binary interactions recommended)
    OUTPATH : Path
        Output path and prefix
    N_LOCI : integer
        Number of loci in final set
    KEEPLIST : Fasta path
        Contains primers that must be included in panel (Default: None)
    SEED : CSV path
        Initial primer set to start with from previous multiplex_wormhole run
    MIN_DIMER : Numeric
        Minimum # dimers to calculate acceptance probabilities for
    MAX_DIMER : Numeric
        Maximum # dimers to calculate acceptance probabilities for
    DECAY_RATE : numeric (0-1)
        Multiplier for negative exponential decay of temperatures (higher values = slower decay)
    T_INIT : Numeric >= 0
        Starting temperature for fixed schedule simulated annealing
    T_FINAL : Numeric < T_INIT
        Ending temperature for filxed schedule simulated annealing
    BURNIN : Numeric
        Iterations used to sample cost space before calculating temperatures
    DIMER_ADJ : Numeric (0-1)
        Proportion of max dimer load to consider when setting T_INIT (high values = more 'bad' changes accepted)
    PROB_ADJ : Numeric (1-Inf)
        Parameters used to adjust dimer acceptance probabilities in exponential model (higher = lower 'bad' changes accepted)
    -------
    Iteratively optimizes a multiplex primer set of given size to minimize
            predicted dimer content. NOTE: This script relies on having a large ratio of 
            loci to choose from relative to the number of loci needed.
    Outputs : CSVs of selected primer pairs + dimer loads
    """
    ## OPTION a: SIMULATED ANNEALING WITH FIXED SCHEDULE
    if T_INIT is not None and T_FINAL is not None:
        pass
    
    ## OPTION B: SIMULATED ANNEALING WITH ADAPTIVE SCHEDULE
    else:
        ## STEP 1: Calculate min and max observed dimer load changes (if not provided)
        if MAX_DIMER is None and MIN_DIMER is None:
            # check for necessary files before proceeding...
            if not [PRIMER_FASTA is not None and DIMER_SUMS is not None and DIMER_TABLE is not None and N_LOCI is not None]:
                raise Exception("InputError! If MAX_DIMER and MIN_DIMER are not specified, paths to PRIMER_FASTA, DIMER_SUMS, "+
                                "DIMER_TABLE and N_LOCI must be provided. Please provide either A) MAX_DIMER and MIN_DIMER or B)"+
                                " PRIMER_FASTA, DIMER_SUMS, DIMER_TABLE, and N_LOCI parameters.")
            else:
                ## Read in IDs and primers
                print("Reading in inputs..........")
                primer_loci, primer_seqs, primer_IDs, primer_pairs = LoadPrimers(PRIMER_FASTA)
            
                # read in dimer info
                #print("Reading in primer dimer counts........")
                dimer_table, dimer_primerIDs, dimer_loci, dimer_tallies, dimer_pairID = LoadDimers(
                    DIMER_SUMS, DIMER_TABLE)
            
                # read in keeplist info
                if KEEPLIST is not None:
                    keeplist_loci, keeplist_pairs = LoadPrimers(KEEPLIST, True)
                    n_keeplist = len(keeplist_pairs)
                else:
                    keeplist_pairs = []
                    keeplist_loci = []
                    n_keeplist = 0
                
                if SEED is None:
                    ## A) Choose initial primer set based on pseudo-Greedy algorithm
                    # Choose initial set of loci based on loci with minimum dimer counts
                    # (Hopefully choosing an initial set this way, rather than randomly, will mean fewer iterations are needed)
                    print("Generating initial primer set........")
                    # make list of unique loci
                    # convert to list because new versions of random.sample won't be able to handle sets...
                    uniq_loci = list(set(primer_loci))
                    nloci = len(uniq_loci)
                    # grab best primer pairs for each locus
                    best_primer_pairs = BestPrimers(uniq_loci, dimer_loci, dimer_tallies, dimer_primerIDs, keeplist_pairs)
                    # if there are fewer loci than desired, use all of them
                    if nloci < N_LOCI:
                        print("WARNING: Fewer loci passed filtering than desired in panel")
                        print("# loci used: " + str(nloci))
                        initial_pairs = best_primer_pairs
                    else:
                        if len(keeplist_pairs) > 0:
                            # remove these options
                            for k in list(best_primer_pairs.keys()):
                                if k in set(keeplist_pairs):
                                    best_primer_pairs.pop(k)
                            # grab N primer pairs with min dimer count (accounting for space filled by keeplist pairs)
                            initial_pairs = dict(sorted(best_primer_pairs.items(), key=lambda x: x[1])[:N_LOCI-n_keeplist])
                            # append all keeplist pairs to initial pairs
                            for pair in set(keeplist_pairs):
                                pair_dimers = [dimer_tallies[x] for x in range(len(dimer_tallies)) if dimer_primerIDs[x] == pair][0]
                                initial_pairs.update({pair: pair_dimers})
                        else:
                            initial_pairs = dict(
                                sorted(best_primer_pairs.items(), key=lambda x: x[1])[:N_LOCI])
                
                ## If a SEED file is provided, use this as the initial primer set....
                else:
                    print("Loading initial primer set from seed file........")
                    initial_pairs = []
                    with open(SEED, 'r') as file:
                        lines = file.readlines()
                        for line in lines:
                            initial_pairs.append(lines[0])
            
                # grab locus IDs for primer pairs
                current_pairIDs = list(initial_pairs.keys())
                current_locusIDs = [GetLocusID(pair) for pair in current_pairIDs]
                # get initial list of allowed alternative primer pairs (i.e., primer pairs for loci not currently in set)
                allowed_loci = [uniq_loci[i] for i in range(
                    len(uniq_loci)) if uniq_loci[i] not in current_locusIDs]
                allowed_idx = list(
                    filter(lambda x: primer_loci[x] in allowed_loci, range(len(primer_loci))))
                allowed_pairs = [primer_pairs[i] for i in allowed_idx]
                allowed_pairs = list(set(allowed_pairs))
                # calculate # primer dimers per pair for current set of primers
                primerset_dimers, nonset_dimers = SubsetDimerTable(
                    current_pairIDs, dimer_table, dimer_pairID, True)
                curr_dimer_totals = CalcTotalDimers(
                    primerset_dimers)  # totals for current primers
                curr_total = sum(curr_dimer_totals.values())
                
                
                ## B) Sample cost space
                print("Sampling dimer load changes.......")
                # A) sample x iterations
                change = []
                i = 0
                while i < BURNIN:
                    # make a new set by randomly swapping a primer pair
                    # newset: 1) replaced ID, 2) new ID, 3) current pair list
                    swap_id, new_id, new_pairIDs = MakeNewSet(current_pairIDs, allowed_pairs, curr_dimer_totals, nonset_dimers, None,
                                                              dimer_primerIDs, dimer_table, dimer_pairID, primer_pairs, primer_loci,
                                                              random=True, keeplist=keeplist_loci)
            
                    # compare newSet to original set
                    comparison, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total = compareSets(new_pairIDs, curr_total, swap_id, new_id, dimer_table, dimer_pairID)
            
                    # if newSet is worse, make note of change value
                    if comparison > 0:
                        change.append(comparison)
                        # repeat
                        i += 1
                MIN_DIMER = min(change)
                MAX_DIMER = max(change)
                print(".....Maximum dimer load observed "+ str(MAX_DIMER))
                print(".....Minimum dimer load observed "+ str(MIN_DIMER))


        ## STEP 2: Set adaptive temperature schedule based on cost calcs
        print("Setting temperature schedule...")
        if T_FINAL is None:
            # set initial and final temps
            T_FINAL = 0
        if T_INIT is None:
            # initial temp should accept most changes
            T_INIT = MIN_DIMER + DIMER_ADJ * (MAX_DIMER)#assuming MIN_DIMER=0
    
    ## PLOT 1: GENERATE TEMPERATURE SCHEDULE ACROSS 100 ITERATIONS
    # Define range of dimer loads to use in calculations
    dimers = list(range(round(MIN_DIMER), round(MAX_DIMER)+1))
    temps = [T_INIT]
    i=0
    while i < 100:
        i+=1
        new_temp = (T_INIT-T_FINAL)*DECAY_RATE**i+T_FINAL
        temps.append(new_temp)
    print(".....Initial temp for adaptive simulated annealing: "+str(T_INIT))
    print(".....Final temp for adaptive simulated annealing: "+str(T_FINAL))
    print("Plotting temperature schedule over 100 iterations")
    plt.figure() # new plotting figure
    plt.plot(temps)
    plt.ylabel("Temperature")
    plt.xlabel("% of Iterations")
    plt.title("Temperature Schedule")
    plt.savefig(OUTPATH+"_TemperatureSchedule.png")
    
    
    ## PLOT 2: DIMER ACCEPTANCE PROBABILITY BY SA TEMPERATURE
    print("Plotting acceptance probabilities for dimers across temperatures...")
    plt.figure() # new plotting figure
    # grab color gradient
    reds = mpl.colormaps['Reds']
    cspace = 1 / len(dimers) # rescale dimer length to 0-1
    # plot each dimer load
    costs = []
    i = 0
    n = 0 # starting value for color scale
    while i<len(dimers):
        d = dimers[i]
        #update color scale value
        n+=cspace
        c = reds(n)
        costs.append([math.exp(-d/x) for x in temps])
        plt.plot(temps, costs[i], color=c, label=str(d)+' dimers')
        i+=1
    # rescale axes and add labels
    plt.title("Dimer Acceptance By Simulated Annealing Temperature")
    plt.xlabel("Temperature decay")
    plt.ylabel("Probability of accepting dimers")
    plt.legend(loc='upper right')
    plt.xlim(T_INIT, T_FINAL)
    #plt.ylim(0, 1)
    plt.savefig(OUTPATH+"_DimerAcceptanceByTemp.png")
    
    
    ## PLOT 3: DIMER ACCEPTANCE PROBABILITY BY % ITERATIONS
    plt.figure()
    i=0; n=0
    while i<len(costs):
        d=dimers[i]
        n+=cspace
        c = reds(n)
        plt.plot(costs[i], color=c, label=str(d)+' dimers')
        i+=1
    plt.title("Dimer Acceptance As Simulated Annealing Progresses")
    plt.xlabel("% of Iterations")
    plt.ylabel("Probability of acceptance")
    plt.legend(loc="upper right")
    plt.savefig(OUTPATH+"_DimerAcceptanceByIteration.png")    



if __name__=="__main__":
    main(sys.arv[1], sys.arv[2], sys.arv[3], sys.argv[4], sys.argv[5], sys.arv[6], sys.argv[7],
         sys.arv[8], sys.arv[9], sys.arv[10], sys.argv[11], sys.argv[12], sys.argv[13], sys.argv[14], sys.argv[15])