#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  2 18:31:45 2025

@author: maggiehallerud
"""

# load dependencies
import sys
import os
import pandas as pd
import math
import random as rand

# load other multiplex wormhole modules
sys.path.append(os.path.dirname(__file__))
import optimize_primers as op
from plot_SA_temps import main as plotSAtemps



def main(PRIMER_FASTA, DIMER_SUMS, DIMER_TABLE, OUTPATH, START_N, KEEPLIST=None, N_RUNS=1, 
         BURNIN=100, ITERATIONS=3000, SIMPLE=3000, PARTITIONS=1000, T_INIT=2.0, T_FINAL=0.1, 
         DECAY_RATE=0.98, DIMER_ADJ=0.1, PROB_ADJ=2, VERBOSE=True, MAKEPLOT=False):
    """

    Parameters
    ----------
    PRIMER_FASTA : TYPE
        DESCRIPTION.
    DIMER_SUMS : TYPE
        DESCRIPTION.
    DIMER_TABLE : TYPE
        DESCRIPTION.
    N_LOCI : TYPE
        DESCRIPTION.
    KEEPLIST : TYPE, optional
        DESCRIPTION. The default is None.
    BURNIN : TYPE, optional
        DESCRIPTION. The default is 100.
    ITERATIONS : TYPE, optional
        DESCRIPTION. The default is 3000.
    SIMPLE : TYPE, optional
        DESCRIPTION. The default is 3000.
    T_INIT : TYPE, optional
        DESCRIPTION. The default is 2.0.
    T_FINAL : TYPE, optional
        DESCRIPTION. The default is 0.1.
    DECAY_RATE : TYPE, optional
        DESCRIPTION. The default is 0.98.
    DIMER_ADJ : TYPE, optional
        DESCRIPTION. The default is 0.1.
    PROB_ADJ : TYPE, optional
        DESCRIPTION. The default is 2.
    VERBOSE : TYPE, optional
        DESCRIPTION. The default is True.

    Raises
    ------
    InputError
        DESCRIPTION.
    Exception
        DESCRIPTION.

    Returns
    -------
    curr_total : TYPE
        DESCRIPTION.

    """
    ## -----------------------STEP 0A: CHECK INPUT PARAMETERS -----------------------##
    op.CheckInputNumber(START_N, "N_LOCI","--")
    op.CheckInputNumber(ITERATIONS, "ITERATIONS", "3000")
    op.CheckInputNumber(SIMPLE, "SIMPLE", "3000")
    
    if ITERATIONS>0:
        op.CheckInputNumber(BURNIN, "BURNIN", "100")
        op.CheckInputNumber(DECAY_RATE, "DECAY_RATE", "0.98")
        op.CheckInputNumber(T_INIT, "T_INIT", "0.1")
        op.CheckInputNumber(T_FINAL, "T_FINAL", "set adaptively, 2 recommended otherwise.")
        op.CheckInputNumber(PARTITIONS,"PARTITIONS", "1000")
        op.CheckInputNumber(DIMER_ADJ, "DIMER_ADJ", "0.1")
        op.CheckInputNumber(PROB_ADJ, "PROB_ADJ", "2")
        
        if T_INIT is None and T_FINAL is None and (BURNIN is None or BURNIN==0):
            raise op.InputError("Either BURNIN or (T_INIT and T_FINAL) need to be provided to run simulated annealing ITERATIONS.")
        if DECAY_RATE<0 or DECAY_RATE>1:
            raise op.InputError("DBECAY_RATE must be between 0 and 1.")
        if DIMER_ADJ<0 or DIMER_ADJ>1:
            raise op.InputError("DIMER_ADJ must be between 0 and 1.")
        if PROB_ADJ<0:
            raise op.InputError("PROB_ADJ must be positive.")
        if (T_INIT is not None and T_FINAL is not None) and T_INIT<T_FINAL:
            raise op.InputError("T_INIT must be greater than T_FINAL.")
        if T_FINAL is not None and T_FINAL<0:
            raise op.InputError("T_FINAL must be a positive number.")    
    
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

    # read in IDs and primers
    primer_loci, primer_seqs, primer_IDs, primer_pairs = op.LoadPrimers(PRIMER_FASTA)

    # read in dimer info
    #print("Reading in primer dimer counts........")
    #dimer_table, dimer_primerIDs, dimer_loci, dimer_tallies, dimer_pairID = LoadDimers(DIMER_SUMS, DIMER_TABLE)
    dimer_table = pd.read_csv(DIMER_TABLE)
    dimer_sums = pd.read_csv(DIMER_SUMS)
    
    # read in keeplist info
    if KEEPLIST is not None:
        keeplist_loci, keeplist_seqs, keeplist_IDs, keeplist_pairs = op.LoadPrimers(KEEPLIST)
        n_keeplist = len(set(keeplist_pairs))
    else:
        keeplist_pairs = []
        keeplist_seqs = []
        keeplist_IDs = []
        n_keeplist = 0
        
    # make list of unique loci
    # convert to list because new versions of random.sample won't be able to handle sets...
    uniq_loci = list(set(primer_loci))
    nloci = len(uniq_loci)
    
    # start at N_LOCI=2 if no starting N provided
    if START_N is None:
        N_LOCI = 2
    else:
        N_LOCI = START_N
    
    # if dimer load already higher than 0, recommend lowering START_N
    initial_pairIDs, allowed, set_dimers, nonset_dimers, curr_dimer_totals, curr_total = SelectInitialPrimers(N_LOCI, uniq_loci, dimer_sums, keeplist_pairs)
    if curr_total>0:
            raise Exception("Dimer load higher than 0 at initial step- rerun with reduced START_N.")
    
    # proceed as long as dimer load is equal to 0
    while curr_total==0:
        print("")
        print("-----------------"+str(N_LOCI)+"-----------------------")
        # try N_RUNS times to meet curr_total==0
        run=1
        while run <= N_RUNS:
            # make initial set
            initial_pairIDs, allowed, set_dimers, nonset_dimers, curr_dimer_totals, curr_total = SelectInitialPrimers(N_LOCI, uniq_loci, dimer_sums, keeplist_pairs)
            # run simulated annealing algorithm
            if ITERATIONS>0:
                print("Starting simulating annealing........")
                if run==1:
                    print("Setting temperature schedule for "+str(N_LOCI)+".........")
                    T_init, T_final, T_space, T_iter = runBurnin(initial_pairIDs, allowed, 
                                                                curr_dimer_totals, curr_total, dimer_table, 
                                                                nonset_dimers, primer_pairs, primer_loci, 
                                                                primer_IDs, primer_seqs, OUTPATH+"_"+str(N_LOCI)+"_", 
                                                                keeplist_IDs, keeplist_seqs, keeplist_pairs)
                curr_pairIDs, allowed, nonset_dimers, curr_dimer_totals, curr_total, costs = \
                    runSA(initial_pairIDs, allowed, curr_dimer_totals, nonset_dimers, 
                          curr_total, dimer_table, primer_pairs, primer_loci, primer_IDs, uniq_loci,
                          primer_seqs, keeplist_seqs, keeplist_pairs, keeplist_IDs, 
                          T_init, T_final, T_iter, T_space, OUTPATH, VERBOSE,
                          ITERATIONS, DECAY_RATE, PROB_ADJ)
                if curr_total==0:
                    N_LOCI +=1
                    continue
            # run greedy algorithm
            if SIMPLE>0:
                print("Starting simple iterative improvement.......")
                curr_pairIDs, curr_dimer_totals, curr_total, costs, nonset_dimers = \
                    runGreedy(curr_pairIDs, keeplist_pairs, allowed, 
                              dimer_table, nonset_dimers, curr_dimer_totals, curr_total, 
                              primer_pairs, primer_loci, primer_IDs,
                              primer_seqs, keeplist_seqs, keeplist_IDs, costs, ITERATIONS,
                              OUTPATH, SIMPLE, VERBOSE)
            # STOP if this is equivalent to the number of loci available
            if nloci < (N_LOCI-n_keeplist):
                raise Exception("Maxed out number of loci available at "+str(N_LOCI))
            # if curr_total>0, try another run
            if curr_total>0:
                run +=1
            else:                
                N_LOCI +=1 # update N_LOCI
                continue # return to start of loop

        ## EXPORT AT EACH STEP ##
        print("EXPORTING OPTIMIZED PRIMER SET....")
        op.ExportCSVs(OUTPATH, curr_dimer_totals, primer_pairs, curr_pairIDs, primer_IDs, 
                      primer_seqs, keeplist_IDs, keeplist_seqs, costs)



def SelectInitialPrimers(N_LOCI, uniq_loci, primer_loci, primer_pairs, dimer_sums, dimer_table, keeplist_pairs):
    # Choose initial set of loci based on loci with minimum dimer counts
    # grab best primer pairs for each locus
    best_primer_pairs = op.BestPrimers(uniq_loci, dimer_sums, keeplist_pairs)
    # if KEEPLIST provided, add additional "best" loci to keeplist to fill out panel
    if len(keeplist_pairs) > 0:
        # remove these options from best for each locus
        for k in list(best_primer_pairs.keys()):
            if k in set(keeplist_pairs):
                best_primer_pairs.pop(k)
        # grab N primer pairs with min dimer count (accounting for space filled by keeplist pairs)
        n_keeplist = len(keeplist_pairs)
        initial_pairs = dict(sorted(best_primer_pairs.items(), key=lambda x: x[1])[:N_LOCI-n_keeplist])
        # grab random subset of pairs
        #initial_keys = rand.sample(best_primer_pairs.items(), N_LOCI-n_keeplist)
        # append all keeplist pairs to initial pairs
        keeplist_dimers = {dimer_sums['Pair1'][x]: dimer_sums['0'][x] for x in range(len(dimer_sums['0'])) if dimer_sums['Pair1'][x] in set(keeplist_pairs)}
        initial_pairs.update(keeplist_dimers)
    # otherwise, with no keeplist, just select the N_LOCI "best" primer pairs to start with
    else:
        initial_pairs = dict(sorted(best_primer_pairs.items(), key=lambda x: x[1])[:N_LOCI])

    # grab locus IDs for primer pairs
    current_pairIDs = list(initial_pairs.keys())
    current_locusIDs = [op.GetLocusID(pair) for pair in current_pairIDs]
    # calculate # primer dimers per pair for current set of primers
    primerset_dimers, nonset_dimers = op.SubsetDimerTable(current_pairIDs, dimer_table, True)
    curr_dimer_totals = op.CalcTotalDimers(primerset_dimers)  
    # totals for current primers
    curr_total = sum(curr_dimer_totals.values())#total

    # get initial list of allowed alternative primer pairs (i.e., primer pairs for loci not currently in set)
    allowed_loci = [uniq_loci[i] for i in range(
        len(uniq_loci)) if uniq_loci[i] not in current_locusIDs]
    allowed_idx = list(
        filter(lambda x: primer_loci[x] in allowed_loci, range(len(primer_loci))))
    allowed_pairs = [primer_pairs[i] for i in allowed_idx]
    allowed_pairs = list(set(allowed_pairs))

    print("     Initial dimer load: "+str(curr_total))
    return current_pairIDs, allowed_pairs, primerset_dimers, nonset_dimers, curr_dimer_totals, curr_total



def runBurnin(initial_pairIDs, allowed_pairs, curr_dimer_totals, curr_total, dimer_table,
              nonset_dimers, primer_pairs, primer_loci, primer_IDs, primer_seqs, OUTPATH,
              keeplist_IDs, keeplist_seqs, keeplist_pairs, MAKEPLOT=False, DECAY_RATE=0.98, PROB_ADJ=2,
              T_INIT=2, T_FINAL=0.1, BURNIN=100, DIMER_ADJ=0.1, ITERATIONS=5000, PARTITIONS=1000):
    # set change log
    change = []
    if (T_INIT is None or T_FINAL is None) and BURNIN>0:
        print("")
        print("SETTING TEMPERATURE SCHEDULE.....")
        print("Running burnin to sample dimer cost space:")
        i = 0
        while i < BURNIN:
            # make a new set by randomly swapping a primer pair
            # newset: 1) replaced ID, 2) new ID, 3) current pair list
            swap_id, new_id, new_pairIDs = op.MakeNewSet(initial_pairIDs, allowed_pairs, curr_dimer_totals, nonset_dimers, [],
                                                         primer_pairs, primer_loci, primer_IDs, primer_seqs, keeplist_IDs, keeplist_seqs, 
                                                         [], [], random=True, keeplist=keeplist_pairs)
            # compare newSet to original set
            comparison, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total = op.compareSets(new_pairIDs, curr_total, swap_id, new_id, dimer_table)
    
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
    
    return T_init, T_final, Tspace, T_iter



def runSA(current_pairIDs, allowed_pairs, curr_dimer_totals, nonset_dimers, 
          curr_total, dimer_table, primer_pairs, primer_loci, primer_IDs, uniq_loci,
          primer_seqs, keeplist_seqs, keeplist_pairs, keeplist_IDs, 
          T_init, T_final, T_iter, T_space, OUTPATH, VERBOSE=False,
          ITERATIONS=5000, DECAY_RATE=0.98, PROB_ADJ=2):
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
            swap_id, new_id, new_pairIDs = op.MakeNewSet(current_pairIDs, allowed_pairs, curr_dimer_totals, nonset_dimers, [],
                                                         primer_pairs, primer_loci,
                                                         OUTPATH, primer_IDs, primer_seqs, keeplist_IDs, keeplist_seqs, costs, [],
                                                         random=True, keeplist=keeplist_pairs, n_iter=n_iter, Temp=Temp, curr_total=curr_total)
            # if e^(-change/temp) < random number (0,1), accept the change
            comparison, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total = op.compareSets(new_pairIDs, curr_total, swap_id, new_id, dimer_table)
            SAvalue = math.exp(-PROB_ADJ*comparison/Temp)
            if SAvalue > rand.uniform(0,1):
                update = op.updateSet(swap_id, new_id, new_pairIDs, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total, verbose=False)
                current_pairIDs, curr_total, curr_dimer_totals, primerset_dimers, nonset_dimers = update #parse output into components
                allowed_pairs = op.UpdateAllowedPairs(swap_id, allowed_pairs, primer_loci, primer_pairs, new_id)
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
            # progress tracking
            if VERBOSE:
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
    current_locusIDs = [op.GetLocusID(pair) for pair in current_pairIDs]
    # get initial list of allowed alternative primer pairs (i.e., primer pairs for loci not currently in set)
    allowed_loci = [uniq_loci[i] for i in range(
        len(uniq_loci)) if uniq_loci[i] not in current_locusIDs]
    allowed_idx = list(
        filter(lambda x: primer_loci[x] in allowed_loci, range(len(primer_loci))))
    allowed_pairs = [primer_pairs[i] for i in allowed_idx]
    allowed_pairs = list(set(allowed_pairs))
    return current_pairIDs, allowed_pairs, nonset_dimers, curr_dimer_totals, curr_total, costs



def runGreedy(current_pairIDs, keeplist_pairs, allowed_pairs, dimer_table, nonset_dimers,
              curr_dimer_totals, curr_total, primer_pairs, primer_loci, primer_IDs,
              primer_seqs, keeplist_seqs, keeplist_IDs, costs, n_iter,
              OUTPATH, SIMPLE=3000, VERBOSE=False):
    i = 0
    if len(keeplist_pairs) > 0:
        blockedlist = keeplist_pairs
    else:
        blockedlist = []
    while i < SIMPLE:
        # progress tracking
        if (i)%1000==0 and VERBOSE:
            print("          finished "+str(i)+" iterations")
        # if current total dimers = 0, STOP because we're optimized already
        if curr_total == 0:
            costs.append([n_iter, "NA", curr_total])
            break

        # create new set by replacing current worst pair with a random pair
        newSet = op.MakeNewSet(current_pairIDs, allowed_pairs, curr_dimer_totals, nonset_dimers, blockedlist,
                               primer_pairs, primer_loci, 
                               OUTPATH, primer_IDs, primer_seqs, keeplist_IDs, keeplist_seqs, costs,
                               [], random=False, keeplist=keeplist_pairs, n_iter=n_iter, Temp="NA", curr_total=curr_total)
        if newSet is None:
            comparison = False
        else:
            swap_id, new_best_id, new_pairIDs = newSet
            # test if new set is better than current set (if fewer overall dimers)
            comparison, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total = op.compareSets(new_pairIDs, curr_total, swap_id, new_best_id, dimer_table)

        # if there are fewer dimers in new set than previous, then update the set
        if comparison < 0:
            update = op.updateSet(swap_id, new_best_id, new_pairIDs, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total, verbose=False)
            current_pairIDs, curr_total, curr_dimer_totals, primerset_dimers, nonset_dimers = update #parse output into components
            op.UpdateAllowedPairs(swap_id, allowed_pairs, primer_loci, primer_pairs)
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
            op.UpdateAllowedPairs(swap_id, allowed_pairs_rmv, primer_loci, primer_pairs)
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
                    newSet = op.MakeNewSet(current_pairIDs, allowed_pairs_rmv, curr_dimer_totals, nonset_dimers, blockedlist,
                                           primer_pairs, primer_loci, 
                                           OUTPATH, primer_IDs, primer_seqs, keeplist_IDs, keeplist_seqs, costs,
                                           blockedlist2, random=False, keeplist=keeplist_pairs, n_iter=n_iter, Temp="NA", curr_total=curr_total)
                    if newSet is None:
                        if VERBOSE:
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
                        if VERBOSE:
                            print("     Infinite loop caused by "+new_best_id+"! Removing from swap options.")
                        blockedlist2.append(new_best_id)
                    # compare new set against current set
                    comparison, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total = op.compareSets(new_pairIDs, curr_total, swap_id, new_best_id, dimer_table)
                    # keep new set if it has fewer dimers
                    if comparison < 0:
                        update = op.updateSet(swap_id, new_best_id, new_pairIDs, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total, verbose=False)
                        current_pairIDs, curr_total, curr_dimer_totals, primerset_dimers, nonset_dimers = op.update #parse output into components
                        allowed_pairs_rmv = op.UpdateAllowedPairs(swap_id, allowed_pairs_rmv, primer_loci, primer_pairs)
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
            costs.append([n_iter, "NA", curr_total])
            if VERBOSE:
                print("No new sets can be made with the available primers!")
            break

        # update iteration #
        i+=1
        # progress tracking
        if VERBOSE:
            print("...Current # dimers: "+str(curr_total))
            if i==SIMPLE:
                costs.append([n_iter, "NA", curr_total])
                print("STOPPED - MAX ITERATIONS REACHED!")

    # report the blockedlisted loci
    if VERBOSE:
        print("The following loci could not be improved with available replacements. If primer loads are high, consider removing these from the input file and rerunning.")
        for b in range(len(blockedlist)):
            print("\t\t" + blockedlist[b])
    
    return current_pairIDs, curr_dimer_totals, curr_total, costs



if __name__=="__main__":
    main()