    #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: OPTIMIZE MULTIPLEX PRIMER SET
Purpose: Iteratively optimizes a multiplex primer set of given size to minimize
        dimer load. Optionally, the user can provide a 'KEEPLIST' of primers
        which must be included in the final multiplex.
Dependencies: pandas (developed w/ v1.4.4)
              matplotlib (developed w/ v3.5.2)

Created on Mon Jul 24 07:59:12 2023
@author: maggiehallerud
"""

# load dependencies
import sys
import importlib
import os
import logging
from datetime import datetime
import csv
import pandas as pd
import numpy as np
import random as rand
import math
import matplotlib.pyplot as plt
import argparse

# import plot simulated annealing temps module
sys.path.append(os.path.dirname(__file__))
plotASAtemps = importlib.import_module("plot_ASA_temps")
from helpers.logging_setup import setup_logging



def main(PRIMER_FASTA, DIMER_SUMS, DIMER_TABLE, OUTPATH, N_LOCI, KEEPLIST=None, deltaG=False, SEED=None, VERBOSE=False,
         SIMPLE=5000, ITERATIONS=1000, CYCLES=10, BURNIN=200, DECAY_RATE=0.95, T_INIT=None, T_FINAL=0.001, 
         PROB_ADJ=2, MAKEPLOT=False, RNG=12345):
    """
    PRIMER_FASTA : Contains primer IDs and sequences [FASTA]
    DIMER_SUMS : Sum of interactions per primer pair (binary interactions recommended) [CSV]
    DIMER_TABLE : Pairwise interaction table of primer pairs (binary interactions recommended) [CSV]
    OUTPATH : Output path and prefix [filepath]
    N_LOCI : Number of primer pairs in desired multipled [integer]
    KEEPLIST : Contains primers that MUST be included in the final solution [FASTA; default=None]
    deltaG : Minimize mean deltaG [True] or count of dimers- requires deltaG dimer tables! [Default: False]
    SEED : Initial primer set to start with from previous multiplex_wormhole run [CSV]
    -------SIMULATED ANNEALING PARAMETERS-----
    SIMPLE : # iterations in simple iterative improvement optimization step [default=5000]
    ITERATIONS : # iterations per simulated annealing cycle [Default: 1000]
    CYCLES : # Simulated annealing cycles to run [Default: 10]
    BURNIN : # iterations to sample dimer cost space [integer; default=1000]
    DECAY_RATE : parameter for temperature decay rate in negative exponential [Default=0.95]
        -closer to 1: least conservative, explores more of cost space but adds more dimers
        -closer to 0: most conservative, does not add dimers but also can't overcome local optima
        -recommended to set lower with more iterations, higher with more:
            # 1000 iterations - 0.98
            # 10000 iterations - 0.90
            # 5000 iterations - 0.95
    T_INIT : Starting temperature for simulated annealing [Default: None, i.e. set adaptively]
        -values closer to T_FINAL will reduce dimer acceptance probabilities
    T_FINAL : Ending temperature for simulated annealing [Default: 0.1]
        -0 is equivalent to simple iterative improvement, higher values allow more costs
    PROB_ADJ : Parameter used to adjust dimer acceptance probabilities [Default: 2]
        -Try increasing to 2 or 3 if too many dimers are being accepted during simulated annealing
    VERBOSE : Print progress? [Default: False]
    RNG : Seed for random numbers.
    -------
    Iteratively optimizes a multiplex primer set of given size to minimize
            predicted dimer content. NOTE: This script relies on having a large ratio of 
            loci to choose from relative to the number of loci needed.
    Outputs : CSVs of selected primer pairs + dimer loads
    """
    ## -----------------------STEP 0A: CHECK INPUT PARAMETERS -----------------------##
    ## first check files exist....
    CheckInputFile(PRIMER_FASTA, "PRIMER_FASTA", required=True)
    CheckInputFile(DIMER_SUMS, "DIMER_SUMS", required=True)
    CheckInputFile(DIMER_TABLE, "DIMER_TABLE", required=True)
    CheckInputFile(KEEPLIST, "KEEPLIST")
    CheckInputFile(SEED, "SEED")    
    CheckInputNumber(N_LOCI, "N_LOCI","--")
    CheckInputNumber(ITERATIONS, "ITERATIONS", "1000")
    CheckInputNumber(SIMPLE, "SIMPLE", "5000")
    # check ASA params
    if ITERATIONS>0:
        CheckInputNumber(BURNIN, "BURNIN", "1000")
        CheckInputNumber(DECAY_RATE, "DECAY_RATE", "0.95")
        CheckInputNumber(T_INIT, "T_INIT", "set adaptively, 2 recommended otherwise.")
        CheckInputNumber(T_FINAL, "T_FINAL", "0.01")
        CheckInputNumber(CYCLES, "CYCLES", "10")
        CheckInputNumber(PROB_ADJ, "PROB_ADJ", "2")
        
        if T_INIT is None and T_FINAL is None and (BURNIN is None or BURNIN==0):
            raise InputError("Either BURNIN or (T_INIT and T_FINAL) need to be provided to run simulated annealing ITERATIONS.")
        if DECAY_RATE<0 or DECAY_RATE>1:
            raise InputError("DBECAY_RATE must be between 0 and 1.")
        if PROB_ADJ<0:
            raise InputError("PROB_ADJ must be positive.")
        if (T_INIT is not None and T_FINAL is not None) and T_INIT<T_FINAL:
            raise InputError("T_INIT must be greater than T_FINAL.")
        if T_FINAL is not None and T_FINAL<0:
            raise InputError("T_FINAL must be a positive number.")
    
    ## -----------------------STEP 0B: READ INPUT FILES -----------------------##
    # initialize logging
    logger = setup_logging(OUTPATH+".log", VERBOSE, NAME=OUTPATH)    
    #print("Logging printed to: ")
    #for h in logger.handlers:
    #    print(type(h), getattr(h, OUTPATH+'.log', None))
    # log start time & inputs
    logger.info("START TIME: %s", datetime.now().strftime('%m/%d/%Y %I:%M:%S %p'))
    logger.info("")
    logger.info("optimize_multiplex inputs: ")
    logger.info("     PRIMER_FASTA: %s", PRIMER_FASTA)
    logger.info("     DIMER_SUMS: %s", DIMER_SUMS)
    logger.info("     DIMER_TABLE: %s", DIMER_TABLE)
    logger.info("     OUTPATH: %s", OUTPATH)
    logger.info("     N_LOCI: %s", N_LOCI)
    logger.info("     KEEPLIST: %s", KEEPLIST)
    logger.info("     deltaG: %s", deltaG)
    logger.info("     SEED: %s", SEED)
    logger.info("     VERBOSE: %s", VERBOSE)
    logger.info("     SIMPLE: %s", SIMPLE)
    logger.info("     ITERATIONS: %s", ITERATIONS)
    logger.info("     CYCLES: %s", CYCLES)
    logger.info("     BURNIN: %s", BURNIN)
    logger.info("     DECAY_RATE: %s", DECAY_RATE)
    logger.info("     T_INIT: %s", T_INIT)
    logger.info("     T_FINAL: %s", T_FINAL)
    logger.info("     PROB_ADJ: %s", PROB_ADJ)
    logger.info("     MAKEPLOT: %s", MAKEPLOT)
    logger.info("     RNG: %s", RNG)
    logger.info("")
    
    logger.info("")
    logger.info("READING INPUTS..........")
    # read in IDs and primers
    primer_loci, primer_seqs, primer_IDs, primer_pairs = LoadPrimers(PRIMER_FASTA)

    # read in dimer info
    #dimer_table, dimer_primerIDs, dimer_loci, dimer_tallies, dimer_pairID = LoadDimers(DIMER_SUMS, DIMER_TABLE)
    dimer_table = pd.read_csv(DIMER_TABLE)
    dimer_sums = pd.read_csv(DIMER_SUMS)
    # convert primer names to string (in case names can be confused for float)
    dimer_table['Pair1'] = [str(x) for x in dimer_table['Pair1']]
    dimer_sums['Pair1'] = [str(x) for x in dimer_sums['Pair1']]
    # check that values match those expected
    neg_sums = dimer_sums['0']<0
    neg_table = dimer_table.iloc[:,1:]<0
    if deltaG and neg_sums.sum()<1:
        raise InputError("No negative values in DIMER_SUMS - are you sure this is the deltaG table and not counts?")
    if deltaG and neg_table.sum().sum()<1:
        raise InputError("No negative values in DIMER_TABLE - are you sure this is the deltaG table and not counts?")
    if not deltaG and neg_sums.sum()>0:
        raise InputError("Negative values in DIMER_SUMS - are you sure this is the count table and not deltaG?")
    if not deltaG and neg_table.sum().sum()>0:
        raise InputError("Negative values in DIMER_TABLE - are you sure this is the count table and not deltaG?")
    # make deltaG values positive s.t. increased costs correspond with worst panels/pairs
    if deltaG:
        dimer_sums['0'] = -dimer_sums['0']
        dimer_table.iloc[:,1:] = -dimer_table.iloc[:,1:]
    
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
    
    # check for missing dimer info
    missing_pairs = [keeplist_pairs[x] for x in range(len(keeplist_pairs)) if keeplist_pairs[x] not in dimer_table['Pair1'].to_numpy()]
    if len(missing_pairs)>0:
        raise InputError("Dimer information is missing for keeplist pairs! Rerun MFEprimer dimer prediction after running AddKeeplist2FASTA")
    primer_pairs, primer_loci, primer_seqs, primer_IDs = CheckPrimersInDimerTables(
        dimer_table, primer_pairs, primer_loci, primer_seqs, primer_IDs, logger)
        
    
    ## -----------------------STEP 1: SELECT INITIAL PRIMER SET -----------------------##
    ## STEP 1: Choose initial primer set based on pseudo-Greedy algorithm
    # Choose initial set of loci based on loci with minimum dimer counts
    # (Hopefully choosing an initial set this way, rather than randomly, will mean fewer iterations are needed)
    logger.info("")
    logger.info("GENERATING INITIAL PRIMER SET........")
    # make list of unique loci
    uniq_loci = list(set(primer_loci))# convert to list because new versions of random.sample won't be able to handle sets...
    # grab best primer pairs for each locus
    best_primer_pairs = BestPrimers(uniq_loci, dimer_sums, keeplist_pairs, deltaG)
    nloci = len(uniq_loci)
    if SEED is None:
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
        logger.info("Read %s primers from SEED file", str(n_seed), SEED)

        # check if keeplist pairs missing from seed pairs- warn if so 
        missed = [x for x in keeplist_pairs if x not in seed_pairIDs]
        if len(missed)>0:
            logger.warning("     WARNING! Some KEEPLIST loci are not in the SEED primer set and will be excluded:")
            i=0
            while i<len(missed):
                logger.info("           "+str(missed[i]))
                i+=1 
        
        # check that length is appropriate
        if n_seed < N_LOCI:
            logger.warning("     WARNING: Adding %s primer pairs to SEED to meet target N_LOCI", str(N_LOCI-n_seed))
            nonseed = {k: v for k, v in best_primer_pairs.items() if k not in seed_pairIDs} #non-SEED primer pairs
            nonseed_pairs = dict(sorted(nonseed.items(), key=lambda x: x[1])[:N_LOCI-n_seed])
            initial_pairs = seed_dict.copy()
            initial_pairs.update(nonseed_pairs)
        elif n_seed > N_LOCI:
            logger.warning("     WARNING: Removing %s worst primer pairs from SEED to meet target N_LOCI", str(n_seed-N_LOCI))
            # calculate dimer load per primer pair in current set
            seed_primerset_dimers, seed_nonset_dimers = SubsetDimerTable(seed_pairIDs, dimer_table, True)
            seed_dimer_totals = CalcTotalDimers(seed_primerset_dimers, deltaG) 
            # keep top N_LOCI pairs with lowest dimer load
            initial_pairs = dict(sorted(seed_dimer_totals.items(), key=lambda x: x[1])[:N_LOCI])
        else:
            # 
            initial_pairs = seed_dict

    # grab locus IDs for primer pairs
    current_pairIDs = list(initial_pairs.keys())
    current_pairIDs = [str(x) for x in current_pairIDs]
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
    curr_dimer_totals = CalcTotalDimers(primerset_dimers, deltaG)  
    # totals for current primers
    if deltaG:
        #mean deltaG of dimers in set (negative s.t. increase = worse dimers, then minimize for min value)
        curr_total = np.array(list(curr_dimer_totals.values())).mean()
        logger.info("     Initial mean deltaG: %s", str(-curr_total))
    else:
        #sum of dimers in count
        curr_total = sum(curr_dimer_totals.values())#total
        logger.info("     Initial dimer load: %s", str(curr_total))
    costs = [["Iterations","ASA_Temp", "TotalDimers"]]

    # if there are fewer loci than desired, raise warning and skip optimization
    if len(initial_pairs) < N_LOCI:
        logger.warning("     WARNING: Fewer templates passed filtering than desired in panel."+
                       "Returning primer pairs with min dimer load for each template.")
        logger.warning("     # loci used: %s", str(len(initial_pairs)))
        # set iterations to zero- no optimization possible
        ITERATIONS = 0
        SIMPLE = 0

    if curr_total==0:
        logger.warning("Initial primer set already has 0 predicted dimers! Skipping optimization process.")
        ITERATIONS=0
        SIMPLE=0
        
    # set initial primer pairs as the current minimum, in case no optimization is possible during SA
    min_pairIDs = current_pairIDs
    min_dimers = curr_total
    min_dimer_totals = curr_dimer_totals
    min_primerset_dimers = primerset_dimers
    min_nonset_dimers = nonset_dimers
    
        
    ## ------------------STEP 2: RUN SIMPLE ITERATIVE IMPROVEMENT OPTIMIZATION----------------##
    # i.e., swap loci only if dimer load is reduced
    # only attempt iterative improvement if no final solution found with simulated annealing
    Temp = 0 #SII = SA w. T=0
    i = 0
    if SIMPLE>0:
        logger.info("")
        logger.info("STARTING SIMPLE ITERATIVE IMPROVEMENT.....")
        if deltaG:
            logger.info("Initial mean deltaG: %s", str(-curr_total))
        else:
            logger.info("Initial # dimers: %s", str(curr_total))
        if len(keeplist_pairs) > 0:
            blockedlist = keeplist_pairs
        else:
            blockedlist = []
        while i < SIMPLE:
            i+=1
            # if current total dimers = 0, STOP because we're optimized already
            if curr_total == 0:
                costs.append([i, Temp, curr_total])
                logger.info("Solution with 0 dimers found!")
                break
    
            # create new set by replacing current worst pair with a random pair
            newSet = MakeNewSet(current_pairIDs, allowed_pairs, curr_dimer_totals, nonset_dimers, blockedlist,
                                primer_pairs, primer_loci, 
                                OUTPATH, primer_IDs, primer_seqs, keeplist_IDs, keeplist_seqs, costs,
                                [], random=False, keeplist=keeplist_pairs, n_iter=i, Temp=Temp, curr_total=curr_total, RNG=12345+i)
            if newSet is None:
                comparison = False
            else:
                swap_id, new_best_id, new_pairIDs = newSet
                # test if new set is better than current set (if fewer overall dimers)
                comparison, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total = \
                    compareSets(new_pairIDs, curr_total, swap_id, new_best_id, dimer_table, dimer_sums, deltaG=deltaG)
                #print(swap_id)
                
            # if there are fewer dimers in new set than previous, then update the set
            if comparison < 0:
                update = updateSet(swap_id, new_best_id, new_pairIDs, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total)
                current_pairIDs, curr_total, curr_dimer_totals, primerset_dimers, nonset_dimers = update #parse output into components
                UpdateAllowedPairs(swap_id, allowed_pairs, primer_loci, primer_pairs)
                costs.append([i, Temp, curr_total])
    
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
                                            blockedlist2, random=False, keeplist=keeplist_pairs, n_iter=i, Temp=Temp, curr_total=curr_total, RNG=12345+i+1)
                        if newSet is None:
                            if VERBOSE:
                                logger.debug("     No new sets can be found for %s ! Removing from swap options.", swap_id)
                            blockedlist.append(swap_id)
                            break
                        else:
                            swap_id, new_best_id, new_pairIDs = newSet
                            break
                        # skip back to top of while if new best is same as original curr_worst
                        if new_best_id == orig_swap or new_best_id == orig_best:
                            #print(new_best_id)
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
                            if VERBOSE:
                                logger.info("     Infinite loop caused by %s ! Removing from swap options.", new_best_id)
                            blockedlist2.append(new_best_id)
                        # compare new set against current set
                        comparison, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total = \
                            compareSets(new_pairIDs, curr_total, swap_id, new_best_id, dimer_table, dimer_sums, deltaG=deltaG)
                        # keep new set if it has fewer dimers
                        if comparison < 0:
                            update = updateSet(swap_id, new_best_id, new_pairIDs, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total)
                            current_pairIDs, curr_total, curr_dimer_totals, primerset_dimers, nonset_dimers = update #parse output into components
                            allowed_pairs_rmv = UpdateAllowedPairs(swap_id, allowed_pairs_rmv, primer_loci, primer_pairs)
                            costs.append([i, Temp, curr_total])
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
                        #print("i")
                        #print("Adding "+curr_worst+" to blockedlist")
                        blockedlist.append(swap_id)
                        continue
    
            # stop iterating if no new sets can be made in the previous while loop
            if newSet is None:
                costs.append([i, Temp, curr_total])
                logger.info("No new sets can be made with the available primers!")
                break
    
            # progress tracking
            if (i)%1000==0:
                logger.info("          finished %s iterations -- dimer load: %s", str(i), str(curr_total))
            # progress tracking
            if i==SIMPLE:
                costs.append([i, Temp, curr_total])
                logger.info("STOPPED - MAX ITERATIONS REACHED!")
    
        # report the blockedlisted loci
        logger.info("The following loci could not be improved with available replacements. If primer loads are high, consider removing these from the input file and rerunning.")
        for b in range(len(blockedlist)):
            logger.info("\t\t %s", blockedlist[b])

        logger.info("")
        logger.info("Dimer load after simple iterative improvement %s", curr_total)
    
            
    ## -----------------------STEP 3: RUN SIMULATED ANNEALING OPTIMIZATION-----------------------##
    if curr_total>0:
        if ITERATIONS>0 and CYCLES>0:
            logger.info("")
            logger.info("STARTING ADAPTIVE SIMULATED ANNEALING OPTIMIZATION......")
            logger.info("Setting up temperature schedule.....")
            ## ------------------STEP 2A: Sample cost space for adaptive temperature schedule----------------##
            if (T_INIT is None) and BURNIN>0 and ITERATIONS>0:
                T_init, MIN_DIMER, MAX_DIMER = setTemps(current_pairIDs, allowed_pairs, curr_dimer_totals, nonset_dimers, 
                                 primer_pairs, primer_loci, OUTPATH, primer_IDs, primer_seqs, keeplist_IDs, 
                                 keeplist_seqs, keeplist_pairs=keeplist_pairs, RNG=RNG, 
                                 curr_total=curr_total, dimer_table=dimer_table, dimer_sums=dimer_sums, 
                                 deltaG=deltaG, logger=logger)
            # use input temps, if provided
            else:
                MIN_DIMER=None; MAX_DIMER=None
                if T_INIT is not None: 
                    logger.info("     Using provided T_INIT=%s", T_INIT)                
                    T_init = T_INIT
                else:
                    logger.info("     Using default T_INIT=3.0")
                    T_init = 3.0
            
            T_final = T_FINAL
            logger.info("     Initial temp: %s", str(round(T_init,2)))
            logger.info("     Final temp: %s", str(round(T_final,2)))
            
            
            ## ------------------STEP 2B: SET AND PLOT TEMPERATURE SCHEDULE----------------##
            # partition iteration space across temp range
            Tspace = ITERATIONS/100 # scaling s.t. temperature schedule is always 0-100
    
            if MIN_DIMER is None: MIN_DIMER=1
            if MAX_DIMER is None: MAX_DIMER=8
            if MAKEPLOT:
                logger.info("     Plotting temperature schedule...")
                plotASAtemps.main(OUTPATH=OUTPATH+"_cycle1",
                                  MIN_DIMER=MIN_DIMER,
                                  MAX_DIMER=MAX_DIMER,
                                  DECAY_RATE=DECAY_RATE,
                                  T_INIT=T_init, 
                                  T_FINAL=T_final,
                                  PROB_ADJ=PROB_ADJ,
                                  BURNIN=0, 
                                  deltaG=deltaG)
        
            ## ------------------STEP 2C: RUN SIMULATED ANNEALING ITERATIONS----------------##
            # C) run simulated annealing optimization
            logger.info("")
            logger.info("STARTING SIMULATED ANNEALING OPTIMIZATION.....")
            # initialize
            asa_costs = []
            Temp = T_init
            cycle = 0
            min_dimers = curr_total
            # start running cycles
            while cycle < CYCLES:
                logger.info("Running cycle %s", cycle+1)
                # set new t_init based on sampling cost space
                if T_INIT is None:
                    if cycle==0:
                        Temp = T_init
                    else:
                        Temp, maxd, mind = setTemps(current_pairIDs, allowed_pairs, curr_dimer_totals, nonset_dimers, 
                                                    primer_pairs, primer_loci, OUTPATH, primer_IDs, primer_seqs, keeplist_IDs, 
                                                    keeplist_seqs, keeplist_pairs, RNG, curr_total, dimer_table, dimer_sums, 
                                                    deltaG, logger, BURNIN=200)
                # always use set T_INIT if provided
                else:
                    Temp = T_INIT
                logger.info("     T_init used: %s", Temp)
                # reset iterations for each cycle
                n_iter = 0
                while n_iter <= ITERATIONS:
                    i+=1
                    n_iter+=1
                    # make a random change
                    newSet = MakeNewSet(current_pairIDs, allowed_pairs, curr_dimer_totals, nonset_dimers, [],
                                        primer_pairs, primer_loci,
                                        OUTPATH, primer_IDs, primer_seqs, keeplist_IDs, keeplist_seqs, costs, [],
                                        random=True, keeplist=keeplist_pairs, n_iter=n_iter, Temp=Temp, curr_total=curr_total, RNG=12345+i)
                    if newSet is not None:
                        swap_id, new_id, new_pairIDs = newSet
                        # log change
                        comparison, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total = \
                            compareSets(new_pairIDs, curr_total, swap_id, new_id, dimer_table, dimer_sums, deltaG=deltaG)
                        asa_costs.append([i, Temp, curr_total])
                        # accept change with probability: e^(-change/temp)
                        SAvalue = math.exp(-PROB_ADJ*comparison/Temp)
                        if SAvalue > rand.uniform(0,1):
                            update = updateSet(swap_id, new_id, new_pairIDs, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total)
                            current_pairIDs, curr_total, curr_dimer_totals, primerset_dimers, nonset_dimers = update #parse output into components
                            allowed_pairs = UpdateAllowedPairs(swap_id, allowed_pairs, primer_loci, primer_pairs, new_id)
                            # log updated cost only when changes are made
                            costs.append([i, Temp, curr_total])
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
                                logger.info("Solution with 0 dimers found!")
                                costs.append([i, Temp, curr_total])
                                break
                    # in simulated annealing, no need to keep track of what's been tested
                    # repeat until # iterations at this temp finishes
                    if VERBOSE:
                        logger.debug(curr_total)
                    # progress tracking
                    if (n_iter)%1000==0:
                        logger.info("     finished %s iterations -- dimer load: %s", str(n_iter), str(curr_total))
    
                    # decrease the temp via exponential decay - this makes the algorithm spend less time in 
                    # 'risky' space and more time in productive optimization space
                    Temp = (T_init-T_final)*DECAY_RATE**(n_iter/Tspace)+T_final
                    #print(str(i)+"   "+str(Temp)+"    "+str(comparison))
                    #print("....."+str(int(n_iter))+" iterations")
                    if curr_total == 0:
                        costs.append([i, Temp, curr_total])
                        logger.info("Solution with 0 dimers found!")
                        break
                # Whenever T_iter is hit, move to the next cycle
                cycle+=1
    
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
            

    ## ------------------STEP 4: EXPORT FINAL PRIMER SETS AND TRACE OF DIMER LOAD----------------##
    logger.info("")
    logger.info("EXPORTING OPTIMIZED PRIMER SET....")
    ExportCSVs(OUTPATH, primer_pairs, current_pairIDs, primer_IDs, 
               primer_seqs, keeplist_IDs, keeplist_seqs, costs, SIMPLE)
    
    logger.info("EXPORTING DIMER TABLES FOR PRIMER PAIRS....")
    final_dimers = SubsetDimerTable(current_pairIDs, dimer_table, return_complement=False)
    # export the current dimers and their totals as CSV
    final_dimers.to_csv(OUTPATH+'_dimers.csv')
    
    # end logging
    logger.info("")
    logger.info("END TIME: %s", datetime.now().strftime('%m/%d/%Y %I:%M:%S %p'))
    logging.shutdown()
    
    return curr_total



def MakeNewSet(pairIDs, allowed, curr_dimer_totals, nonset_dimers, blockedlist, 
               primer_pairs, primer_loci, OUTPATH, primer_IDs, primer_seqs, keeplist_IDs, keeplist_seqs, costs,
               blockedlist2=[], random=False, keeplist=[], n_iter=None, Temp=None, curr_total=None, RNG=12345):
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
        ExportCSVs(OUTPATH, primer_pairs, pairIDs, primer_IDs, 
                   primer_seqs, keeplist_IDs, keeplist_seqs, costs, 0)
        raise OptimizationWarning("STOPPED. Keeplist primer pairs are the only ones with dimer loads, therefore further optimization is impossible. (outputs saved).")

    ## IDENTIFY PRIMER PAIR TO BE REPLACED
    # if swapping randomly, choose a random pairID (not including keeplist) to replace
    if random:
        # remove keeplist and blocklist loci from options
        options = [x for x in pairIDs if x not in keeplist]
        options = [x for x in options if x not in blockedlist]
        rand.seed(RNG)
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
            rand.seed(RNG)
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
                rand.seed(RNG)
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
            rand.seed(RNG)
            new_id = rand.choice(allowed_pairs_edit)
        else:
            return None

    ## UPDATE THE PRIMER SET
    update_pairIDs.append(new_id)
    update_pairIDs.remove(swap)
    return swap, new_id, update_pairIDs



def setTemps(current_pairIDs, allowed_pairs, curr_dimer_totals, nonset_dimers, 
         primer_pairs, primer_loci, OUTPATH, primer_IDs, primer_seqs, keeplist_IDs, 
         keeplist_seqs, keeplist_pairs, RNG, curr_total, dimer_table, dimer_sums, 
         deltaG, logger, BURNIN=200):
    logger.info("     Sampling cost space to set initial temp....")
    change = []
    while len(change) < BURNIN:
        i = RNG+len(change)
        # make a new set by randomly swapping a primer pair
        # newset: 1) replaced ID, 2) new ID, 3) current pair list
        newSet = MakeNewSet(current_pairIDs, allowed_pairs, curr_dimer_totals, nonset_dimers, [],
                                                  primer_pairs, primer_loci, OUTPATH, primer_IDs, primer_seqs, keeplist_IDs, 
                                                  keeplist_seqs, [], [], random=True, keeplist=keeplist_pairs, RNG=i)
        if newSet is not None:
            swap_id, new_id, new_pairIDs = newSet
            # compare newSet to original set
            comparison, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total = \
                compareSets(new_pairIDs, curr_total, swap_id, new_id, dimer_table, dimer_sums, deltaG=deltaG)
            change.append(comparison)
        else:
            break
    
    # set temperatures based on burnin results
    neg_changes = [x for x in change if x>0]
    if len(neg_changes)==0: neg_changes = [0.01]
    MEAN_DIMER = np.mean(neg_changes)
    MAX_DIMER = max(neg_changes)
    MIN_DIMER = min(neg_changes)#almost always 1!
    logger.info("     Changes that increase cost: %s", 
                str(round(len(neg_changes)/len(change),3)*100)+"%")
    logger.info("     Mean +cost change: %s", round(MEAN_DIMER,2))
    logger.info("     Max +cost observed: %s", round(MAX_DIMER,2))
    logger.info("     Min +cost observed: %s", round(MIN_DIMER,2))
    # T_init should depend on p(accepting mistake) [based on mean mistake] & 
    # how bad that 'mistake' is
    if deltaG:
        T_INIT = (-2*math.log10(MEAN_DIMER/MAX_DIMER)+1)/10
        if T_INIT<0.03: T_INIT=0.03
    else:
        T_INIT = -2*math.log10(MEAN_DIMER/MAX_DIMER)+2
        if T_INIT < 0.3: T_INIT = 0.3
    
    return [T_INIT, MIN_DIMER, MAX_DIMER]



def compareSets(new_pairIDs, curr_total, swap, new_best_id, dimer_table, dimer_sums, deltaG):
    # calculate # dimers in this new set
    new_primerset_dimers, new_nonset_dimers = SubsetDimerTable(new_pairIDs, dimer_table, True)
    new_dimer_totals = CalcTotalDimers(new_primerset_dimers, deltaG)
    if deltaG:
        # take mean for deltaG comparisons
        new_total = np.array(list(new_dimer_totals.values())).mean()
    else:
        new_total = sum(new_dimer_totals.values())
    comparison = new_total - curr_total  # difference between new set and old
    # if deltaG, use coarse scaling (mean*#)
    if deltaG:
        comparison = comparison*len(new_pairIDs)

    return comparison, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total


def updateSet(swap_id, new_id, new_pairIDs, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total):
    # update all starting values to match this iteration's values
    current_pairIDs = new_pairIDs.copy()
    primerset_dimers, nonset_dimers = new_primerset_dimers, new_nonset_dimers
    curr_dimer_totals = new_dimer_totals
    curr_total = new_total
    return current_pairIDs, curr_total, curr_dimer_totals, primerset_dimers, nonset_dimers


def UpdateAllowedPairs(swap_pair, inlist, primer_loci, primer_pairs, new_pair = None):
    # add alternative primer pairs for swap locus back to list
    swap_pair = str(swap_pair)
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
    if len(swap_pairs)>0:
        for p in swap_pairs:
            inlist.append(p)
    # remove primer pairs for new pair ID
    if new_pair is not None:
        new_pair = str(new_pair)
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


def BestPrimers(loci, dimer_sums, keeplist, deltaG, RNG=12345):
    """loci: """
    """dimer_tally_dict: """
    # grab dimer data
    dimer_tallies = dimer_sums['0']
    dimer_primerIDs = dimer_sums['Pair1']
    dimer_loci = [str(x).split(".")[0] for x in dimer_primerIDs]
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
            # extract primer pair with worst dimer value
            min_idx = list(filter(lambda x: tallies[x] == min(tallies), range(len(tallies))))
            # if there's more than one 'best' option, randomly choose one...
            if len(min_idx) > 1:
                rand.seed(RNG)
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


def CalcTotalDimers(dimerDF, deltaG):
    pairs = list(dimerDF.columns)
    if 'Pair1' in pairs:
        pairs.remove('Pair1')
    outDict = dict()
    for pair in pairs: 
        if deltaG:
            Ndimers = dimerDF[pair].mean()
        else:
            Ndimers = sum(dimerDF[pair])
        outDict.update({pair: Ndimers})
    return outDict



def CheckPrimersInDimerTables(dimer_table, primer_pairs, primer_loci, primer_seqs, primer_IDs, logger):
    missing_pairs = [primer_pairs[x] for x in range(len(primer_pairs)) if primer_pairs[x] not in dimer_table['Pair1'].to_numpy()]
    if len(missing_pairs)>0:
        logger.warning("%s CANDIDATE PRIMERS WERE NOT PRESENT IN DIMER TABLES! These may have been removed if they would duplicate amplification of the keeplist.", len(missing_pairs))
        logger.warning("Candidates with missing dimer info were removed:")
        primer_loci = [primer_loci[x] for x in range(len(primer_loci)) if primer_pairs[x] not in missing_pairs]
        primer_seqs = [primer_seqs[x] for x in range(len(primer_seqs)) if primer_pairs[x] not in missing_pairs]
        primer_IDs = [primer_IDs[x] for x in range(len(primer_IDs)) if primer_pairs[x] not in missing_pairs]
        primer_pairs = [primer_pairs[x] for x in range(len(primer_pairs)) if primer_pairs[x] not in missing_pairs]
    for p in missing_pairs:
        logger.warning("\t%s", str(p))    
    return primer_pairs, primer_loci, primer_seqs, primer_IDs



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
                primerid = line.rstrip()
                primerid = primerid.replace('>', '')
                locus = str(primerid.split(".")[0])
                pair = primerid.replace(".FWD", "").replace(".REV", "")
                primer_loci.append(locus)
                primer_IDs.append(primerid)
                primer_pairs.append(pair)
            else:
                line = line.rstrip()
                primer_seqs.append(line)
    if keeplist:
        primer_pairs = list(set(primer_pairs))
        primer_loci = []
        for i in primer_pairs:
            locus = str(i.split(".")[0])
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
            dimer_primerIDs.append(str(linesplit[0]))
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



def ExportCSVs(OUTPATH, primer_pairs, current_pairIDs, primer_IDs, primer_seqs, 
               keeplist_IDs, keeplist_seqs, costs, ITERATIONS):
    # extract info for selected primer pairs
    current_pairs_index = list(
        filter(lambda x: primer_pairs[x] in current_pairIDs, range(len(primer_pairs))))
    outpairs = [primer_IDs[x] for x in current_pairs_index]
    outseqs = [primer_seqs[x] for x in current_pairs_index]
    # add keeplist loci
    if len(keeplist_IDs)>0:
        outpairs = outpairs + keeplist_IDs
        outseqs = outseqs + keeplist_seqs    
    # export selected multiplex to CSV
    with open(OUTPATH+'_primers.csv', 'w') as file:
        file.write("PrimerID,Sequence\n")
        for i in range(len(outpairs)):
            try:
                file.write(outpairs[i]+","+outseqs[i]+"\n")
            except Exception:
                pass
    
    if len(costs)>2:
        ## Plot cost trace across iterations
        iterations = [float(costs[i][0]) for i in range(1,len(costs))]
        dimercosts = [float(costs[i][2]) for i in range(1,len(costs))]
        plt.figure()
        plt.vlines(x=int(ITERATIONS), ymin=0, ymax=max(dimercosts),linestyles="dashed",colors="black")
        plt.plot(iterations, dimercosts)
        plt.title("Trace of Dimer Load")
        plt.xlabel("Iterations")
        plt.ylabel("Total Dimer Load")
        plt.savefig(OUTPATH+"_DimerLoad.png")
        
        ## export cost changes
        with open(OUTPATH+'_costsTrace.csv', 'w') as file:
            for line in costs:
                file.write(str(line[0])+","+str(line[1])+","+str(line[2])+"\n")



def CheckInputFile(file, name, required=False):
    if required or file is not None:
        if not os.path.exists(file):
            raise InputError(name+" file could not be found!")



class InputError(Exception):
    pass

class OptimizationWarning(Exception):
    pass



def parse_args():
    # initialize argparser
    parser = argparse.ArgumentParser()
    # add required arguments
    parser.add_argument("-f", "--primer_fasta", type=str, required=True)
    parser.add_argument("-d", "--dimer_sums", type=str, required=True)
    parser.add_argument("-t", "--dimer_table", type=str, required=True)
    parser.add_argument("-o", "--outpath", type=str, required=True)
    parser.add_argument("-n", "--nloci", type=int, required=True)
    # add optional arguments
    parser.add_argument("-k", "--keeplist", type=str, default=None)
    parser.add_argument("-e", "--seed", type=str, default=None)
    parser.add_argument("-s", "--simple", type=int, default=5000)
    parser.add_argument("-i", "--iter", type=int, default=1000)
    parser.add_argument("-c", "--cycles", type=int, default=10)
    parser.add_argument("-b", "--burnin", type=int, default=200)
    parser.add_argument("-r", "--decay_rate", type=float, default=0.95)
    parser.add_argument("-x", "--temp_init", type=float, default=None)
    parser.add_argument("-l", "--temp_final", type=float, default=0.1)
    parser.add_argument("-a", "--prob_adj", type=float, default=2)
    # add flags
    parser.add_argument("-g", "--deltaG", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-m", "--makeplot", action="store_true")
    return parser.parse_args()



if __name__=="__main__":
    # parse command line arguments
    args = parse_args()
    # run main
    main(PRIMER_FASTA = args.primer_fasta, 
         DIMER_SUMS = args.dimer_sums, 
         DIMER_TABLE = args.dimer_table, 
         OUTPATH = args.outpath, 
         N_LOCI = args.nloci, 
         KEEPLIST = args.keeplist, 
         deltaG = args.deltaG, 
         SEED = args.seed, 
         VERBOSE = args.verbose,
         SIMPLE = args.simple, 
         ITERATIONS = args.iter, 
         CYCLES = args.cycles,
         BURNIN = args.burnin, 
         DECAY_RATE = args.decay_rate, 
         T_INIT = args.temp_init, 
         T_FINAL = args.temp_final, 
         PROB_ADJ = args.prob_adj, 
         MAKEPLOT = args.makeplot)