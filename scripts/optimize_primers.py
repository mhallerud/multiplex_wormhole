#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: OPTIMIZE MULTIPLEX PRIMER SET
Created on Mon Jul 24 07:59:12 2023
@author: maggiehallerud

Dependencies: sys, random (should come preinstalled with Python)

Purpose: Iteratively optimizes a multiplex primer set of given size to minimize
        predicted dimer content. NOTE: This script relies on having a large ratio of 
        loci to choose from relative to the number of loci needed.

Inputs: Fasta file of primer inputs, sums of primer pair interactions and pairwise table
        of primer pair interactions.
Outputs: CSV of selected primer pairs and CSV with expected dimer loads per pair.

"""
# objective function = sumcost + missingLoci*penalty
# penalty = N_LOCI / 10 (i.e., a primer can only NOT be added if it forms dimers with >10% of the existing set)

# STEP 1: Greedy algorithm (ish)
# select the units with lowest overall cost

# STEP 2: Adaptive simulated annealing
# A) sample x iterations & set initial and final temps
# final temp (i.e., target temp) = minChange
# initial temp = minChange + 0.1 * (maxChange - minChange)
# min and max based on least and most 'bad' changes made (i.e., lowest and highest increase in dimer load observed)
# initial temp should accept almost any change
# B) run simulated annealing
# make a random change
# if e^(-change/temp) < random number (0,1), accept the change
# decrease the temp with each accepted change (10,000 decreases generally recommended in marxan)
# repeat until final temp is reached
# in simulated annealing, no need to keep track of what's been tested

# STEP 3: Use iterative improvement via random swapping to find the local optima
# only test the same PU once




# load dependencies
import os
import sys
import csv
import pandas as pd
import random as rand
import math
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(__file__))
from plot_SA_temps import main as plotSAtemps



def main(PRIMER_FASTA, DIMER_SUMS, DIMER_TABLE, OUTPATH, N_LOCI, KEEPLIST=None, SEED=None, VERBOSE=False,
         SIMPLE=1000, ITERATIONS=2000, BURNIN=100, DECAY_RATE=0.98, T_INIT=None, T_FINAL=None, PARTITIONS=1000, 
         DIMER_ADJ=0.1, PROB_ADJ=2, MAKEPLOT=True):
    """
    PRIMER_FASTA : Contains primer IDs and sequences [FASTA]
    DIMER_SUMS : Sum of interactions per primer pair (binary interactions recommended) [CSV]
    DIMER_TABLE : Pairwise interaction table of primer pairs (binary interactions recommended) [CSV]
    OUTPATH : Output path and prefix [filepath]
    N_LOCI : Number of primer pairs in desired multipled [integer]
    KEEPLIST : Contains primers that MUST be included in the final solution [FASTA; default=None]
    SEED : Initial primer set to start with from previous multiplex_wormhole run (CSV)
    -------SIMULATED ANNEALING PARAMETERS-----
    SIMPLE : # iterations in simple iterative improvement optimization step [default=5000]
    ITERATIONS : # iterations in simulated annealing optimization step [default=10000]
        -increase iterations for complex problems
        -set to 0 to skip simulated annealing step
    BURNIN : # iterations to sample dimer cost space [integer; default=100]
    DECAY_RATE : parameter for temperature decay rate in negative exponential [proportion; default=0.98]
        -closer to 1: least conservative, explores more of cost space but adds more dimers
        -closer to 0: most conservative, does not add dimers but also can't overcome local optima
        -recommended to set lower with more iterations, higher with more:
            # 1000 iterations - 0.98
            # 10000 iterations - 0.90
            # 5000 iterations - 0.95
    T_INIT : Starting temperature for simulated annealing (default: None, i.e. set adaptively)
        -values closer to T_FINAL will reduce dimer acceptance probabilities
    T_FINAL : Ending temperature for simulated annealing (default: 0.1)
        -0 is equivalent to simple iterative improvement, higher values allow more costs
    PARTITIONS : # partitions for temperature space [default=1000]
        -fewer partitions means algorithm stays in each temperature space longer
        -overriden if iterations < partitions
    DIMER_ADJ : proportion of dimer load to consider when setting SA temperatures [proportion; default=0.1]
        -values closer to 1 will allow many dimers to be accepted
        -values closer to 0 will reject dimers but fail to overcome local optima
    PROB_ADJ : Parameter used to adjust dimer acceptance probabilities (default: 1)
        -Try increasing to 2 or 3 if too many dimers are being accepted during simulated annealing
    -------
    Iteratively optimizes a multiplex primer set of given size to minimize
            predicted dimer content. NOTE: This script relies on having a large ratio of 
            loci to choose from relative to the number of loci needed.
    Outputs : CSVs of selected primer pairs + dimer loads
    """
    
    ## -----------------------STEP 0A: CHECK INPUT PARAMETERS -----------------------##
    CheckInputNumber(N_LOCI, "N_LOCI","--")
    CheckInputNumber(ITERATIONS, "ITERATIONS", "2000")
    CheckInputNumber(SIMPLE, "SIMPLE", "1000")
    
    if ITERATIONS>0:
        CheckInputNumber(BURNIN, "BURNIN", "100")
        CheckInputNumber(DECAY_RATE, "DECAY_RATE", "0.98")
        CheckInputNumber(T_INIT, "T_INIT", "0.1")
        CheckInputNumber(T_FINAL, "T_FINAL", "set adaptively, 2 recommended otherwise.")
        CheckInputNumber(PARTITIONS,"PARTITIONS", "1000")
        CheckInputNumber(DIMER_ADJ, "DIMER_ADJ", "0.1")
        CheckInputNumber(PROB_ADJ, "PROB_ADJ", "2")
        
        if T_INIT is None and T_FINAL is None and (BURNIN is None or BURNIN==0):
            raise InputError("Either BURNIN or (T_INIT and T_FINAL) need to be provided to run simulated annealing ITERATIONS.")
        if DECAY_RATE<0 or DECAY_RATE>1:
            raise InputError("DBECAY_RATE must be between 0 and 1.")
        if DIMER_ADJ<0 or DIMER_ADJ>1:
            raise InputError("DIMER_ADJ must be between 0 and 1.")
        if PROB_ADJ<0:
            raise InputError("PROB_ADJ must be positive.")
        if (T_INIT is not None and T_FINAL is not None) and T_INIT<T_FINAL:
            raise InputError("T_INIT must be greater than T_FINAL.")
        if T_FINAL is not None and T_FINAL<0:
            raise InputError("T_FINAL must be a positive number.")
    
    if SEED is not None and SEED==OUTPATH+"_primers.csv":
        raise InputError("OUTPATH_primers.csv will overwrite SEED file- rename one.")
    
    
    ## -----------------------STEP 0B: READ INPUT FILES -----------------------##
    # check for files first
    print("READING INPUTS..........")
    if not os.path.exists(PRIMER_FASTA):
        raise Exception("CHECK INPUTS! PRIMER_FASTA <"+PRIMER_FASTA+"> file not found.")
    if not os.path.exists(DIMER_SUMS):
        raise Exception("CHECK INPUTS! DIMER_SUMS <"+DIMER_SUMS+"> file not found.")
    if not os.path.exists(DIMER_TABLE):
        raise Exception("CHECK INPUTS! DIMER_TABLE <"+DIMER_TABLE+"> file not found.")
    if KEEPLIST is not None and not os.path.exists(KEEPLIST):
        raise Exception("CHECK INPUTS! KEEPLIST <"+KEEPLIST+"> file not found.")
    if SEED is not None and not os.path.exists(SEED):
        raise Exception("CHECK INPUTS! SEED <"+SEED+"> file not found.")

    # read in IDs and primers
    primer_loci, primer_seqs, primer_IDs, primer_pairs = LoadPrimers(
        PRIMER_FASTA)

    # read in dimer info
    #print("Reading in primer dimer counts........")
    #dimer_table, dimer_primerIDs, dimer_loci, dimer_tallies, dimer_pairID = LoadDimers(DIMER_SUMS, DIMER_TABLE)
    dimer_table = pd.read_csv(DIMER_TABLE)
    dimer_sums = pd.read_csv(DIMER_SUMS)
    
    # read in keeplist info
    if KEEPLIST is not None:
        keeplist_loci, keeplist_seqs, keeplist_IDs, keeplist_pairs = LoadPrimers(KEEPLIST)
        n_keeplist = len(set(keeplist_pairs))
    else:
        keeplist_pairs = []
        keeplist_loci = []
        keeplist_seqs = []
        keeplist_IDs = []
        n_keeplist = 0
    

    ## -----------------------STEP 1: SELECT INITIAL PRIMER SET -----------------------##
    ## STEP 1: Choose initial primer set based on pseudo-Greedy algorithm
    # Choose initial set of loci based on loci with minimum dimer counts
    # (Hopefully choosing an initial set this way, rather than randomly, will mean fewer iterations are needed)
    print("")
    print("GENERATING INITIAL PRIMER SET........")
    # make list of unique loci
    # convert to list because new versions of random.sample won't be able to handle sets...
    uniq_loci = list(set(primer_loci))
    nloci = len(uniq_loci)
    # grab best primer pairs for each locus
    best_primer_pairs = BestPrimers(uniq_loci, dimer_sums, keeplist_pairs)
    if SEED is None:
        # if there are fewer loci than desired, use all of them
        if nloci < (N_LOCI-n_keeplist):
            print("     WARNING: Fewer loci passed filtering than desired in panel")
            print("     # loci used: " + str(nloci))
            initial_pairs = best_primer_pairs
        else:
            # if KEEPLIST provided, add additional "best" loci to keeplist to fill out panel
            if len(keeplist_pairs) > 0:
                # remove these options from best for each locus
                for k in list(best_primer_pairs.keys()):
                    if k in set(keeplist_pairs):
                        best_primer_pairs.pop(k)
                # grab N primer pairs with min dimer count (accounting for space filled by keeplist pairs)
                initial_pairs = dict(sorted(best_primer_pairs.items(), key=lambda x: x[1])[:N_LOCI-n_keeplist])
                # grab random subset of pairs
                #initial_keys = rand.sample(best_primer_pairs.items(), N_LOCI-n_keeplist)
                # append all keeplist pairs to initial pairs
                keeplist_dimers = {dimer_sums['Pair1'][x]: dimer_sums['0'][x] for x in range(len(dimer_sums['0'])) if dimer_sums['Pair1'][x] in set(keeplist_pairs)}
                initial_pairs.update(keeplist_dimers)
            # otherwise, with no keeplist, just select the N_LOCI "best" primer pairs to start with
            else:
                initial_pairs = dict(sorted(best_primer_pairs.items(), key=lambda x: x[1])[:N_LOCI])
    
    ## If a SEED file is provided, use this as the initial primer set....
    else:
        # grab primer pair IDs from SEED
        seed_pairIDs = []
        with open(SEED, 'r') as file:
            lines = file.readlines()
            for line in lines[1:]: #skip header
                line = line.split(",")
                seed_pairIDs.append(line[0].replace(".FWD","").replace(".REV",""))
        # convert to dictionary with total dimer load
        seed_pairIDs = list(set(seed_pairIDs)) #grab unique values since there will be 2 records per primer pair
        seed_dict = {dimer_sums['Pair1'][x]: dimer_sums['0'][x] for x in range(len(dimer_sums['0'])) if dimer_sums['Pair1'][x] in seed_pairIDs}
        n_seed = len(seed_pairIDs)
        print("Read "+str(n_seed)+" primers from SEED file <"+SEED+">")

        # check if keeplist pairs missing from seed pairs- warn if so 
        missed = [x for x in keeplist_pairs if x not in seed_pairIDs]
        if len(missed)>0:
            print("     WARNING! Some KEEPLIST loci are not in the SEED primer set and will be excluded:")
            i=0
            while i<len(missed):
                print("           "+str(missed[i]))
                i+=1 
        
        # check that length is appropriate
        if n_seed < N_LOCI:
            print("     Adding "+str(N_LOCI-n_seed)+" primer pairs to SEED to meet target N_LOCI")
            nonseed = {k: v for k, v in best_primer_pairs.items() if k not in seed_pairIDs} #non-SEED primer pairs
            nonseed_pairs = dict(sorted(nonseed.items(), key=lambda x: x[1])[:N_LOCI-n_seed])
            initial_pairs = seed_dict.copy()
            initial_pairs.update(nonseed_pairs)
        elif n_seed > N_LOCI:
            print("     Removing "+str(n_seed-N_LOCI)+" worst primer pairs from SEED to meet target N_LOCI")
            # calculate dimer load per primer pair in current set
            seed_primerset_dimers, seed_nonset_dimers = SubsetDimerTable(seed_pairIDs, dimer_table, True)
            seed_dimer_totals = CalcTotalDimers(seed_primerset_dimers) 
            # keep top N_LOCI pairs with lowest dimer load
            initial_pairs = dict(sorted(seed_dimer_totals.items(), key=lambda x: x[1])[:N_LOCI])
        else:
            # 
            initial_pairs = seed_dict

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
    primerset_dimers, nonset_dimers = SubsetDimerTable(current_pairIDs, dimer_table, True)
    curr_dimer_totals = CalcTotalDimers(primerset_dimers)  
    # totals for current primers
    curr_total = sum(curr_dimer_totals.values())#total
    print("     Initial dimer load: "+str(curr_total))
    if curr_total==0:
        print("Initial primer set already has 0 predicted dimers! Skipping optimization process.")
        BURNIN=0
        ITERATIONS=0
        SIMPLE=0
        
    # set initial primer pairs as the current minimum, in case no optimization is possible during SA
    min_pairIDs = current_pairIDs
    min_dimers = curr_total
    min_dimer_totals = curr_dimer_totals
    min_primerset_dimers = primerset_dimers
    min_nonset_dimers = nonset_dimers

    
    ## -----------------------STEP 2: RUN SIMULATED ANNEALING OPTIMIZATION-----------------------##
    ## ------------------STEP 2A: Sample cost space for adaptive temperature schedule----------------##
    if ITERATIONS>0:
        change = []
        if (T_INIT is None or T_FINAL is None) and BURNIN>0 and ITERATIONS>0:
            print("")
            print("SETTING TEMPERATURE SCHEDULE.....")
            print("Running burnin to sample dimer cost space:")
            i = 0
            while i < BURNIN:
                # make a new set by randomly swapping a primer pair
                # newset: 1) replaced ID, 2) new ID, 3) current pair list
                swap_id, new_id, new_pairIDs = MakeNewSet(current_pairIDs, allowed_pairs, curr_dimer_totals, nonset_dimers, [],
                                                          primer_pairs, primer_loci, 
                                                          OUTPATH, primer_IDs, primer_seqs, keeplist_IDs, keeplist_seqs, [], [],
                                                          random=True, keeplist=keeplist_pairs)
                # compare newSet to original set
                comparison, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total = compareSets(new_pairIDs, curr_total, swap_id, new_id, dimer_table)
        
                # if newSet is worse, make note of change value
                if comparison > 0:
                    change.append(comparison)
                    # repeat
                    i += 1
            print("     Max change observed: "+str(max(change)))    
            # set temperatures based on burnin results
            if T_FINAL is None:
                T_final = 0.1 #float(min(change)) # default setting
            if T_INIT is None:
                T_init = min(change) + DIMER_ADJ * (max(change) - min(change)) # adaptive simulated annealing
        # use input temps, if provided
        else:
            if T_INIT is not None: 
                print("Using provided T_INIT")                
                T_init = T_INIT
            else:
                print("Using default T_INIT")
                T_init = 2.0
            if T_FINAL is not None:
                print("Using provided T_FINAL")                
                T_final = T_FINAL
            else:
                print("Using default T_FINAL")
                T_final=0.1
        print("     Initial temp: "+str(T_init))
        print("     Final temp: "+str(T_final))
        
        ## ------------------STEP 2B: SET AND PLOT TEMPERATURE SCHEDULE----------------##
        # partition iteration space across temp range
        if ITERATIONS > PARTITIONS:
            Tspace = (T_init - T_final)/PARTITIONS
            T_iter = round(Tspace*ITERATIONS)
        else:
            print("    NOTE: Fewer partitions provided than iterations, +\
                  so the temperature schedule is set to change at each iteration.")
            Tspace = (T_init - T_final)/ITERATIONS
            T_iter = 1    
        if len(change)==0:
            MIN_DIMER=1
            MAX_DIMER=10
        else:
            MIN_DIMER=min(change)
            MAX_DIMER=max(change)
        if MAKEPLOT:
            print("     Plotting temperature schedule...")
            plotSAtemps(OUTPATH=OUTPATH,
                        MIN_DIMER=MIN_DIMER,
                        MAX_DIMER=MAX_DIMER,
                        DECAY_RATE=DECAY_RATE,
                        T_INIT=T_init, 
                        T_FINAL=T_final, 
                        DIMER_ADJ=DIMER_ADJ,
                        PROB_ADJ=PROB_ADJ,
                        BURNIN=0)
    
        ## ------------------STEP 2C: RUN SIMULATED ANNEALING ITERATIONS----------------##
        # C) run simulated annealing optimization
        print("")
        print("STARTING SIMULATED ANNEALING OPTIMIZATION.....")
        # initialize
        Temp = T_init
        n_iter = 0
        step = 0
        min_dimers = curr_total
        # start tracking costs
        costs = [["Iterations","SA_Temp", "TotalDimers"]]
        costs.append([n_iter, Temp, curr_total]) #log for cost values
        # start iterating
        while n_iter < ITERATIONS:
            i = 1  # reset iterations for each temperature
            while i <= T_iter:
                # make a random change
                swap_id, new_id, new_pairIDs = MakeNewSet(current_pairIDs, allowed_pairs, curr_dimer_totals, nonset_dimers, [],
                                                          primer_pairs, primer_loci,
                                                          OUTPATH, primer_IDs, primer_seqs, keeplist_IDs, keeplist_seqs, costs, [],
                                                          random=True, keeplist=keeplist_pairs, n_iter=n_iter, Temp=Temp, curr_total=curr_total)
                # if e^(-change/temp) < random number (0,1), accept the change
                comparison, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total = compareSets(new_pairIDs, curr_total, swap_id, new_id, dimer_table)
                SAvalue = math.exp(-PROB_ADJ*comparison/Temp)
                if SAvalue > rand.uniform(0,1):
                    update = updateSet(swap_id, new_id, new_pairIDs, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total, verbose=False)
                    current_pairIDs, curr_total, curr_dimer_totals, primerset_dimers, nonset_dimers = update #parse output into components
                    allowed_pairs = UpdateAllowedPairs(swap_id, allowed_pairs, primer_loci, primer_pairs, new_id)
                    # log updated cost only when changes are made
                    costs.append([n_iter, Temp, curr_total])
                    # keep track of set with minimum dimers
                    # this allows the best set to be kept while the algorithm continues to explore more of the space
                    if curr_total <= min_dimers:
                        min_pairIDs = current_pairIDs
                        min_dimers = curr_total
                        min_dimer_totals = curr_dimer_totals
                        min_primerset_dimers = primerset_dimers 
                        min_nonset_dimers = nonset_dimers
                    # STOP if 0 dimers in set
                    if curr_total == 0:
                        costs.append([n_iter, Temp, curr_total])
                        break
                # in simulated annealing, no need to keep track of what's been tested
                # repeat until # iterations at this temp finishes
                i+=1
                if VERBOSE:
                    print(curr_total)
                # progress tracking
                if (n_iter+i)%1000==0:
                    print("          finished "+str(n_iter+i)+" iterations -- dimer load: "+str(curr_total))

            # decrease the temp via exponential decay - this makes the algorithm spend less time in 
            # 'risky' space and more time in productive optimization space
            step+=1
            Temp = (T_init-T_final)*DECAY_RATE**step+T_final
            n_iter+=T_iter
            #print("....."+str(int(n_iter))+" iterations")
            if curr_total == 0:
                costs.append([n_iter, Temp, curr_total])
                print("Solution with 0 dimers found!")
                break
        
        # proceed with best set found during simulated annealing- which isn't necessarily the final set
        current_pairIDs = min_pairIDs
        curr_total = min_dimers
        curr_dimer_totals = min_dimer_totals
        primerset_dimers = min_primerset_dimers 
        nonset_dimers = min_nonset_dimers
        # grab locus IDs for current primer pairs
        current_pairIDs = list(curr_dimer_totals.keys())
        current_locusIDs = [GetLocusID(pair) for pair in current_pairIDs]
        # get initial list of allowed alternative primer pairs (i.e., primer pairs for loci not currently in set)
        allowed_loci = [uniq_loci[i] for i in range(
            len(uniq_loci)) if uniq_loci[i] not in current_locusIDs]
        allowed_idx = list(
            filter(lambda x: primer_loci[x] in allowed_loci, range(len(primer_loci))))
        allowed_pairs = [primer_pairs[i] for i in allowed_idx]
        allowed_pairs = list(set(allowed_pairs))
        
        ## ------------------STEP 2D: SAVE SA OUTPUTS AS CHECKPOINT----------------##
        # export the current dimers and their totals as CSV
        ExportCSVs(OUTPATH+"_SA", curr_dimer_totals, primer_pairs, current_pairIDs, 
                   primer_IDs, primer_seqs, keeplist_IDs, keeplist_seqs, costs)
    

    ## ------------------STEP 3: RUN SIMPLE ITERATIVE IMPROVEMENT OPTIMIZATION----------------##
    # i.e., swap loci only if dimer load is reduced
    # only attempt iterative improvement if no final solution found with simulated annealing
    if SIMPLE>0:
        print("STARTING GREEDY OPTIMIZATION.....")
        print("Initial # dimers: " + str(curr_total))
        i = 0
        if len(keeplist_pairs) > 0:
            blockedlist = keeplist_pairs
        else:
            blockedlist = []
        while i < SIMPLE:
            # progress tracking
            if (i)%1000==0:
                print("          finished "+str(i)+" iterations")
            # if current total dimers = 0, STOP because we're optimized already
            if curr_total == 0:
                costs.append([n_iter, Temp, curr_total])
                print("Solution with 0 dimers found!")
                break
    
            # create new set by replacing current worst pair with a random pair
            newSet = MakeNewSet(current_pairIDs, allowed_pairs, curr_dimer_totals, nonset_dimers, blockedlist,
                                primer_pairs, primer_loci, 
                                OUTPATH, primer_IDs, primer_seqs, keeplist_IDs, keeplist_seqs, costs,
                                [], random=False, keeplist=keeplist_pairs, n_iter=n_iter, Temp=Temp, curr_total=curr_total)
            if newSet is None:
                comparison = False
            else:
                swap_id, new_best_id, new_pairIDs = newSet
                # test if new set is better than current set (if fewer overall dimers)
                comparison, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total = compareSets(new_pairIDs, curr_total, swap_id, new_best_id, dimer_table)
    
            # if there are fewer dimers in new set than previous, then update the set
            if comparison < 0:
                update = updateSet(swap_id, new_best_id, new_pairIDs, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total, verbose=False)
                current_pairIDs, curr_total, curr_dimer_totals, primerset_dimers, nonset_dimers = update #parse output into components
                UpdateAllowedPairs(swap_id, allowed_pairs, primer_loci, primer_pairs)
                costs.append([n_iter+i, 0, curr_total])
    
            # otherwise, loop through remaining replacement options
            # this avoids retesting the same options for this replacement
            else:
                # remove tested pair from allowed lists
                allowed_pairs_rmv = allowed_pairs.copy()
                try:
                    allowed_pairs_rmv.remove(new_best_id)
                except:
                    pass
                # add indices for the other primer pairs for the current worst locus to the allowed list
                UpdateAllowedPairs(swap_id, allowed_pairs_rmv, primer_loci, primer_pairs)
                # avoid trying the same values...
                orig_swap = swap_id
                orig_best = new_best_id
                # loop through, allowing all other replacement options to be tested
                x = len(allowed_pairs_rmv)
                blockedlist2 = []
                while x > 0:
                    if x > 0:
                        # set these to avoid infinite loops
                        prev_worst = swap_id
                        prev_best = new_best_id
                        # make new set
                        newSet = MakeNewSet(current_pairIDs, allowed_pairs_rmv, curr_dimer_totals, nonset_dimers, blockedlist,
                                            primer_pairs, primer_loci, 
                                            OUTPATH, primer_IDs, primer_seqs, keeplist_IDs, keeplist_seqs, costs,
                                            blockedlist2, random=False, keeplist=keeplist_pairs, n_iter=n_iter, Temp=Temp, curr_total=curr_total)
                        if newSet is None:
                            print("     No new sets can be found for "+swap_id+"! Removing from swap options.")
                            blockedlist.append(swap_id)
                            break
                        else:
                            swap_id, new_best_id, new_pairIDs = newSet
                        # skip back to top of while if new best is same as original curr_worst
                        if new_best_id == orig_swap or new_best_id == orig_best:
                            if x > 1:
                                try:
                                    allowed_pairs_rmv.remove(new_best_id)
                                except:
                                    pass
                                x = len(allowed_pairs_rmv)
                                continue  # go back to while
                            else:
                                break
                        # add to blockedlist2 if stuck in an infinite loop
                        if new_best_id == prev_best and swap_id == prev_worst:
                            print("     Infinite loop caused by "+new_best_id+"! Removing from swap options.")
                            blockedlist2.append(new_best_id)
                        # compare new set against current set
                        comparison, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total = compareSets(new_pairIDs, curr_total, swap_id, new_best_id, dimer_table)
                        # keep new set if it has fewer dimers
                        if comparison < 0:
                            update = updateSet(swap_id, new_best_id, new_pairIDs, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total, verbose=False)
                            current_pairIDs, curr_total, curr_dimer_totals, primerset_dimers, nonset_dimers = update #parse output into components
                            allowed_pairs_rmv = UpdateAllowedPairs(swap_id, allowed_pairs_rmv, primer_loci, primer_pairs)
                            break  # exit this loop if a replacement was found
                        else:
                            try:
                                allowed_pairs_rmv.remove(new_best_id)
                            except:
                                pass
                            try:
                                allowed_pairs_rmv.remove(swap_id)
                            except:
                                pass
                    # otherwise, loop through remaining replacement options
                    # this avoids retesting the same options for this replacement
                    x = len(allowed_pairs_rmv)  # update
                    # when options run out for replacing this pair, add it to blockedlist to avoid trying to replace
                    if x <= 1:
                        #print("Adding "+curr_worst+" to blockedlist")
                        blockedlist.append(swap_id)
    
            # stop iterating if no new sets can be made in the previous while loop
            if newSet is None:
                costs.append([n_iter, Temp, curr_total])
                print("No new sets can be made with the available primers!")
                break
    
            # update iteration #
            i+=1
            # progress tracking
            if VERBOSE:
                print("...Current # dimers: "+str(curr_total))
            if i==SIMPLE:
                costs.append([n_iter, Temp, curr_total])
                print("STOPPED - MAX ITERATIONS REACHED!")
    
        # report the blockedlisted loci
        print("The following loci could not be improved with available replacements. If primer loads are high, consider removing these from the input file and rerunning.")
        for b in range(len(blockedlist)):
            print("\t\t" + blockedlist[b])
    

        ## ------------------STEP 4: EXPORT FINAL PRIMER SETS AND TRACE OF DIMER LOAD----------------##
        print("")
        print("EXPORTING OPTIMIZED PRIMER SET....")
        ExportCSVs(OUTPATH, curr_dimer_totals, primer_pairs, current_pairIDs, primer_IDs, 
                   primer_seqs, keeplist_IDs, keeplist_seqs, costs)
                    
    return curr_total



def MakeNewSet(pairIDs, allowed, curr_dimer_totals, nonset_dimers, blockedlist, 
               primer_pairs, primer_loci, OUTPATH, primer_IDs, primer_seqs, keeplist_IDs, keeplist_seqs, costs,
               blockedlist2=[], random=False, keeplist=[], n_iter=None, Temp=None, curr_total=None):
    # create copy of original list, otherwise links them and affects both lists when using update/remove
    update_pairIDs = pairIDs.copy()
    allowed_pairs_edit = allowed.copy()
    nonset_dimers_edit = nonset_dimers.copy()
    
    ## CHECK WHETHER FURTHER OPTIMIZATION IS POSSIBLE
    dimersub = {k: v for k, v in curr_dimer_totals.items() if k not in keeplist}
    worst = sorted(dimersub.items(), key=lambda x: x[1])
    # raise message if all primer pairs outside of the keeplist already have 0 dimers
    worst_dimers = [x[1] for x in worst]
    worst_sum = sum(worst_dimers)
    if worst_sum==0:
        costs.append([n_iter, Temp, curr_total])
        ExportCSVs(OUTPATH, curr_dimer_totals, primer_pairs, pairIDs, primer_IDs, 
                   primer_seqs, keeplist_IDs, keeplist_seqs, costs)
        raise OptimizationWarning("STOPPED. Keeplist primer pairs are the only ones with dimer loads, therefore further optimization is impossible. (outputs saved).")

    ## IDENTIFY PRIMER PAIR TO BE REPLACED
    # if swapping randomly, choose a random pairID (not including keeplist) to replace
    if random:
        # remove keeplist and blocklist loci from options
        options = [x for x in pairIDs if x not in keeplist]
        options = [x for x in options if x not in blockedlist]
        swap = rand.choice(options)  
    # otherwise, replace the current worst (most dimers) pair in the set
    else:
        if len(dimersub)>0:
            worst_pair = worst[-1]
            swap = worst_pair[0]
        else:
            return None
        
        # if the chosen ID is in the keeplist or the blockedlist of loci that can't be fixed, 
        # then choose the next max value
        while swap in blockedlist:
            worst.remove(worst_pair)
            if len(worst)>0:    
                worst_pair = worst[-1]
                swap = worst_pair[0]
            else:
                return None # this means no options could be found

    # add other primers for the swap locus back to the allowed list
    allowed_pairs_edit = UpdateAllowedPairs(swap, allowed_pairs_edit, primer_loci, primer_pairs)
    # Remove any primer pairs that cause infinite loops
    for pair in blockedlist2:
        try:
            allowed_pairs_edit.remove(pair)
        except:
            pass

    ## CHOOSE REPLACEMENT PRIMER PAIR
    if len(allowed_pairs_edit) > 0:
        if random:
            new_id = rand.choice(allowed_pairs_edit)
        else:
            ## ID primer pair outside of set that would have the smallest dimer load in the context of the current set
            ## The main purpose of this is to speed up the optimization (as opposed to just choosing a random pair, many)
            ## of which may not be very good in any context). This method also doesn't allow any flexibility if the min option doesn't work.
            # remove swap ID from dimer table
            nonset_dimers_edit.loc[~nonset_dimers_edit['Pair1'].isin([swap])]
            #nonset_dimers_edit.pop(swap)
            # subset dimer table with rows = current primer pairs, columns = candidates
            candidate_dimers = nonset_dimers.loc[:,nonset_dimers.columns.isin(allowed_pairs_edit)]
            #nonset_pairIDs = list(set(dimer_pairID) - set(pairIDs))
            #allowed_indx = [x for x in range(len(nonset_pairIDs)) if nonset_pairIDs[x] in allowed_pairs_edit]
            #candidate_ids =[x for x in nonset_pairIDs if x in allowed_pairs_edit]
            #candidate_dimers = dict()
            #for pair in nonset_dimers_edit.keys():
            #    pairDict = nonset_dimers_edit[pair]
            #    subPair = [pairDict[i] for i in allowed_indx]
            #    candidate_dimers.update({pair: subPair})
            # grab column sums (total dimer load for nonset/candidate dimers)
            candidate_totals = candidate_dimers.sum(axis=0)
                #candidate_totals = list(map(sum, zip(*candidate_dimers.values()))) #get sums across same column of values
            if len(candidate_totals)>0:
                # grab the new pair ID (randomly choose if multiple options)
                cand_min = min(candidate_totals)
                new_best_options = candidate_totals[candidate_totals==cand_min]
                #new_best_options = [candidate_ids[x] for x in range(len(candidate_ids)) if candidate_totals[x]==cand_min]
                new_id = rand.choice(new_best_options.keys())
            else:
                return None
    else:
        return None
    
    # double check that chosen pair isn't already in set... if it is, remove and try again
    while new_id in update_pairIDs:
        if new_id in allowed_pairs_edit:
            allowed_pairs_edit.remove(new_id)
        if len(allowed_pairs_edit) > 0:
            new_id = rand.choice(allowed_pairs_edit)
        else:
            return None

    ## UPDATE THE PRIMER SET
    update_pairIDs.append(new_id)
    update_pairIDs.remove(swap)
    return swap, new_id, update_pairIDs


def compareSets(new_pairIDs, curr_total, swap, new_best_id, dimer_table):
    # calculate # dimers in this new set
    new_primerset_dimers, new_nonset_dimers = SubsetDimerTable(new_pairIDs, dimer_table, True)
    new_dimer_totals = CalcTotalDimers(new_primerset_dimers)
    new_total = sum(new_dimer_totals.values())
    comparison = new_total - curr_total  # difference between new set and old
    return comparison, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total


def updateSet(swap_id, new_id, new_pairIDs, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total, verbose=True):
        if verbose:
            print(".........."+swap_id + " replaced with " + new_id)
    # if this primer set has fewer dimers overall, keep it and proceed
    #if new_total < curr_total:
        # TODO: if a primer pair has been tried, remove from current 'allowed' set (so that we don't go in circles...)
        # maybe not necessary though because it all depends on what's still in the set...
        # update all starting values to match this iteration's values
        current_pairIDs = new_pairIDs.copy()
        primerset_dimers, nonset_dimers = new_primerset_dimers, new_nonset_dimers
        curr_dimer_totals = new_dimer_totals
        curr_total = new_total
        return current_pairIDs, curr_total, curr_dimer_totals, primerset_dimers, nonset_dimers


def UpdateAllowedPairs(swap_pair, inlist, primer_loci, primer_pairs, new_pair = None):
    # add alternative primer pairs for swap locus back to list
    swap_locus = swap_pair.split(".")[0]
    swap_idx = list(
        filter(lambda x: primer_loci[x] == swap_locus, range(len(primer_loci))))
    swap_pairs = [primer_pairs[i] for i in swap_idx]
    swap_pairs = (set(swap_pairs))  # extract unique values only
    # make sure current worse isn't in here...
    try:
        swap_pairs.remove(swap_pair)
    except Exception:
        pass
    for p in swap_pairs:
        inlist.append(p)
    # remove primer pairs for new pair ID
    if new_pair is not None:
        new_locus = new_pair.split(".")[0]
        new_idx = list(filter(lambda x: primer_loci[x] == new_locus, range(len(primer_loci))))
        new_pairs = [primer_pairs[i] for i in new_idx]
        new_pairs = (set(new_pairs))  # extract unique values only
        for p in new_pairs:
            try:
                inlist.remove(p)
            except Exception:
                pass
    # return new allowed list
    return inlist


def condition(x, list):
    if x in list:
        return x


def GetLocusID(pairID):
    pairIDsplit = pairID.split(".")
    LocusID = pairIDsplit[0]
    return LocusID


def BestPrimers(loci, dimer_sums, keeplist):
    """loci: """
    """dimer_tally_dict: """
    # grab dimer data
    dimer_tallies = dimer_sums['0']
    dimer_primerIDs = dimer_sums['Pair1']
    dimer_loci = [x.split(".")[0] for x in dimer_primerIDs]
    # replace "best" with keeplist pair if in keeplist
    keeplist_loci = [keeplist[x].split('.')[0] for x in range(len(keeplist))]
    best_primer_pairs = dict()
    for locus in loci:
        if locus in keeplist_loci:
            # assign keeplist pairs as best
            keeplist_indx = list(
                filter(lambda x: keeplist_loci[x] == locus, range(len(keeplist_loci))))
            bestPair = [keeplist[x] for x in keeplist_indx][0]
            # extract dimer tally for this pair
            dimer_idx = list(
                filter(lambda x: dimer_primerIDs[x] == bestPair, range(len(dimer_primerIDs))))
            if len(dimer_idx) == 0:
                min_dimers = float('nan')  # in case pair is missing from DB
            else:
                min_dimers = [dimer_tallies[x] for x in dimer_idx][0]
        else:
            # extract dimer data for this locus
            locus_idx = list(
                filter(lambda x: dimer_loci[x] == locus, range(len(dimer_loci))))
            ids = [dimer_primerIDs[i] for i in locus_idx]
            tallies = [dimer_tallies[i] for i in locus_idx]
            # extract primer pair with min dimer value
            min_idx = list(filter(lambda x: tallies[x] == min(tallies), range(len(tallies))))
            # if there's more than one 'best' option, randomly choose one...
            if len(min_idx) > 1:
                min_idx = [rand.choice(min_idx)]
            # extract the ID for the selected pair
            bestPair = [ids[i] for i in min_idx][0]
            min_dimers = min(tallies)
        # set up dictionary with locus: (primer pair id, min # dimers)
        best_primer_pairs.update({bestPair: min_dimers})
    return best_primer_pairs


def SubsetDimerTable(pairs, dimer_table, return_complement=False):
    """primer_pairs"""
    """dimer_tally_dict"""
    # grab rows, columns for primer pairs in dimer table
    rows = dimer_table['Pair1'].isin(pairs)
    cols = dimer_table.columns.isin(pairs)
    pair1 = dimer_table.columns.get_loc("Pair1")
    cols[pair1] = True # reset to True to include pair field
    # subset dimer table for primer pairs
    sub_dimer_table = dimer_table.loc[rows, cols]
    # subset dimer table - all other primer pairs, with rows as those in current set
    cols[pair1] = False
    not_dimer_table = dimer_table.loc[rows, ~cols] #~ means "not"
    #sub_dimer_table = dict()  # blank dictionary to store output
    #not_dimer_table = dict()
    # grab dictionary indices for the input primer pairs (skip first ID because that's the name column)
    #primer_idx = list(filter(lambda x: dimer_ids[x] in primer_pairs, range(1, len(dimer_ids))))
    # left-adjust by 1 due to empty first column in dimer_ids
    #primer_idx = [i-1 for i in primer_idx]
    #for pair in primer_pairs:
    #    # grab all primer pair interactions between the specified pairs
    #    pair_dict = dimer_table[pair]
    #    sub_pair_dict = [int(pair_dict[i]) for i in primer_idx]
    #    sub_dimer_table.update({pair: sub_pair_dict})
    #    # grab all primer pair interactions between the pair and any pairs not specified
    #    not_pair_dict = [int(pair_dict[i])
    #                     for i in range(len(pair_dict)) if i not in primer_idx]
    #    not_dimer_table.update({pair: not_pair_dict})
    # return subset table (and table with non-subset values)
    if not return_complement:
        return sub_dimer_table
    else:
        return sub_dimer_table, not_dimer_table


def CalcTotalDimers(dimerDF):
    pairs = list(dimerDF.columns)
    if 'Pair1' in pairs:
        pairs.remove('Pair1')
    outDict = dict()
    for pair in pairs: 
        Ndimers = sum(dimerDF[pair])
        outDict.update({pair: Ndimers})
    return outDict


def LoadPrimers(PRIMER_FASTA, keeplist=False):
    # read in files as lists
    primer_loci = []
    primer_IDs = []
    primer_seqs = []
    primer_pairs = []
    # parse sequences from headers
    with open(PRIMER_FASTA, "r") as file:
        lines = file.readlines()
        for line in lines:
            if ">" in line:
                line = line.rstrip()
                line = line.replace('>', '')
                locus = line.split(".")[0]
                pair = line.replace(".FWD", "").replace(
                    ".FW", "").replace(".REV", "")
                primer_loci.append(locus)
                primer_IDs.append(line)
                primer_pairs.append(pair)
            else:
                line = line.rstrip()
                primer_seqs.append(line)
    if keeplist:
        primer_pairs = list(set(primer_pairs))
        primer_loci = []
        for i in primer_pairs:
            locus = i.split(".")[0]
            primer_loci.append(locus)
        return primer_loci, primer_pairs
    else:
        return primer_loci, primer_seqs, primer_IDs, primer_pairs


def LoadDimers(DIMER_SUMS, DIMER_TABLE):
    dimer_primerIDs = []
    dimer_loci = []
    dimer_tallies = []
    with open(DIMER_SUMS, 'r') as file:
        lines = file.readlines()
        for line in lines[1:]:  # skip first line because R puts unnecessary header
            # read in line and reformat
            line = line.rstrip()
            linesplit = line.split(",")
            # remove extra quotes from string
            linesplit[0] = linesplit[0].replace('"', '')
            linesplitid = linesplit[0].split(".")
            # extract info
            dimer_primerIDs.append(linesplit[0])
            dimer_loci.append(linesplitid[0])
            dimer_tallies.append(int(linesplit[1]))
    # we'll store this in a dictionary so that it's really easy to subset based on the primer pair ID
    dimer_table = dict()
    dimer_pairID = []
    with open(DIMER_TABLE, 'r', newline="\n") as file:
        reader = csv.reader(file, delimiter=',')
        next(reader)  # first line is header so we're skipping it
        for row in reader:
            dimer_table.update({row[0]: row[1:]})
            dimer_pairID.append(row[0])
    return dimer_table, dimer_primerIDs, dimer_loci, dimer_tallies, dimer_pairID



def CheckInputNumber(X, NAME, DEFAULT):
    if X is not None:
        try:
            xint = int(X)
        except ValueError:
            try:
                xfl = float(X)
            except ValueError:
                raise Exception("CHECK INPUTS! "+NAME+" is not a number. [DEFAULT="+DEFAULT+"]")


def ExportCSVs(OUTPATH, curr_dimer_totals, primer_pairs, current_pairIDs, primer_IDs, primer_seqs, 
               keeplist_IDs, keeplist_seqs, costs):
    # export the current dimers and their totals as CSV
    with open(OUTPATH+'_dimers.csv', 'w') as file:
        file.write("PrimerPairID,DimerLoad\n")
        for key in curr_dimer_totals.keys():
            # line = curr_dimer_totals[]
            #line_str = str(line)[1:-1].replace("'","")
            file.write(key+","+str(curr_dimer_totals[key])+"\n")

    # export selected primers to CSV
    current_pairs_index = list(
        filter(lambda x: primer_pairs[x] in current_pairIDs, range(len(primer_pairs))))
    outpairs = [primer_IDs[x] for x in current_pairs_index]
    outseqs = [primer_seqs[x] for x in current_pairs_index]
    # add keeplist loci
    if len(keeplist_IDs)>0:
        outpairs = outpairs + keeplist_IDs
        outseqs = outseqs + keeplist_seqs    
    with open(OUTPATH+'_primers.csv', 'w') as file:
        file.write("PrimerID,Sequence\n")
        for i in range(len(outpairs)):
            try:
                file.write(outpairs[i]+","+outseqs[i]+"\n")
            except Exception:
                pass
    
    ## Plot cost trace across iterations
    iterations = [int(costs[i][0]) for i in range(1,len(costs))]
    dimercosts = [int(costs[i][2]) for i in range(1,len(costs))]
    plt.figure()
    #plt.vlines(x=ITERATIONS, ymin=0, ymax=max(dimercosts),linestyles="dashed",colors="black")
    plt.plot(iterations, dimercosts)
    plt.title("Trace of Dimer Load")
    plt.xlabel("Iterations")
    plt.ylabel("Total Dimers")
    plt.savefig(OUTPATH+"_DimerLoad.png")

    
    # export cost changes
    with open(OUTPATH+'_ASA_costs.csv', 'w') as file:
        for line in costs:
            file.write(str(line[0])+","+str(line[1])+","+str(line[2])+"\n")



class InputError(Exception):
    pass

class OptimizationWarning(Exception):
    pass


# run main function when called from the command line
if __name__ == "__main__":
    main(sys.argv[1], #PRIMER_FASTA
         sys.argv[2], #DIMER_SUMS
         sys.argv[3], #DIMER_TABLE
         sys.argv[4], #OUTPATH
         int(sys.argv[5]), #N_LOCI
         sys.argv[6], #KEEPLIST
         sys.argv[7], #SEED
         sys.arv[8], #VERBOSE
         sys.argv[9], #SIMPLE
         sys.argv[10], #ITERATIONS
         sys.argv[11], #BURNIN
         sys.argv[12], #DECAY_RATE
         sys.argv[13], #T_INIT
         sys.argv[14], #T_FINAL
         sys.argv[15], #PARTITIONS
         sys.argv[16], #DIMER_ADJ
         sys.argv[17]) #PROB_ADJ
