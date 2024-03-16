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


## SIMULATED ANNEALING PARAMETERS:
# BURNIN (# iterations to sample space)
BURNIN = 100
# ADJUSTMENT (proportion of dimer load rand to consider when setting initial temp)
# - values closer to 1 will accept nearly any change - extremely long run time and very 'bad' changes accepted
# - values closer to 0 will only accept changes very similar to the 'least bad' options - quicker but won't explore much of the space
ADJUSTMENT = 0.1
# PARTITIONs (# partitions for temperature space)
PARTITIONS = 1000
# - fewer partitions means algorithm stays in each temperature space longer
# - partitions are overriden if iterations < partitions
# ITERATIONS = # iterations in adaptive simulated annealing (after temps are set)
ITERATIONS = 5000
# - more iterations: slower but better for optimization
# - fewer iterations may be OK for simpler problems where local minima are not expected to be a problem
# SIMPLE = # iterations in simple iterative improvement
SIMPLE = 500
# DECAY_RATE
DECAY_RATE = 0.95 #must be between 0 and 1 (non-inclusive)
# - closer to 1: algorithm will spend more time at high temperatures, and take more risks (adding more dimers)
# - closer to 0.5: algorithm will spend little time at high temperatures and take few risks
# - closer to 0: algorithm will spend virtually no time at high temperatures
# set lower with more iterations, higher with more:
# 1000 iterations - 0.98
# 10000 iterations - 0.90
# 5000 iterations - 0.95


# load dependencies
#import os
import sys
import csv
import random as rand
import math



def main(PRIMER_FASTA, DIMER_SUMS, DIMER_TABLE, OUTPATH, N_LOCI, WHITELIST=None, SEED=None):
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
    WHITELIST : Fasta path
        Contains whitelist primer IDs and sequences (Default: None)
    SEED : CSV path
        Initial primer set to start with from previous multiplex_wormhole run
    -------
    Iteratively optimizes a multiplex primer set of given size to minimize
            predicted dimer content. NOTE: This script relies on having a large ratio of 
            loci to choose from relative to the number of loci needed.
    Outputs : CSVs of selected primer pairs + dimer loads
    """
    ## STEP 0: Read in IDs and primers
    print("Reading in inputs..........")
    primer_loci, primer_seqs, primer_IDs, primer_pairs = LoadPrimers(
        PRIMER_FASTA)

    # read in dimer info
    #print("Reading in primer dimer counts........")
    dimer_table, dimer_primerIDs, dimer_loci, dimer_tallies, dimer_pairID = LoadDimers(
        DIMER_SUMS, DIMER_TABLE)

    # read in whitelist info
    if WHITELIST is not None:
        whitelist_loci, whitelist_pairs = LoadPrimers(WHITELIST, True)
        n_whitelist = len(whitelist_pairs)
    else:
        whitelist_pairs = []
        n_whitelist = 0
    
    if SEED is None:
        ## STEP 1: Choose initial primer set based on pseudo-Greedy algorithm
        # Choose initial set of loci based on loci with minimum dimer counts
        # (Hopefully choosing an initial set this way, rather than randomly, will mean fewer iterations are needed)
        print("Generating initial primer set........")
        # make list of unique loci
        # convert to list because new versions of random.sample won't be able to handle sets...
        uniq_loci = list(set(primer_loci))
        nloci = len(uniq_loci)
        # grab best primer pairs for each locus
        best_primer_pairs = BestPrimers(
            uniq_loci, dimer_loci, dimer_tallies, dimer_primerIDs, whitelist_pairs)
        # if there are fewer loci than desired, use all of them
        if nloci < N_LOCI:
            print("WARNING: Fewer loci passed filtering than desired in panel")
            print("# loci used: " + str(nloci))
            initial_pairs = best_primer_pairs
        else:
            if len(whitelist_pairs) > 0:
                # remove these options
                for k in list(best_primer_pairs.keys()):
                    if k in set(whitelist_pairs):
                        best_primer_pairs.pop(k)
                # grab N primer pairs with min dimer count (accounting for space filled by whitelist pairs)
                initial_pairs = dict(sorted(best_primer_pairs.items(), key=lambda x: x[1])[:N_LOCI-n_whitelist])
                # append all whitelist pairs to initial pairs
                for pair in set(whitelist_pairs):
                    pair_dimers = [dimer_tallies[x] for x in range(len(dimer_tallies)) if dimer_primerIDs[x] == pair][0]
                    initial_pairs.update({pair: pair_dimers})
            else:
                initial_pairs = dict(
                    sorted(best_primer_pairs.items(), key=lambda x: x[1])[:N_LOCI])
    
    ## If a SEED file is provided, use this as the initial primer set....
    else:
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
    
    
    # STEP 2: Adaptive simulated annealing
    print("Sampling dimer load changes.......")
    # A) sample x iterations
    change = []
    i = 0
    while i < BURNIN:
        # make a new set by randomly swapping a primer pair
        # newset: 1) replaced ID, 2) new ID, 3) current pair list
        swap_id, new_id, new_pairIDs = MakeNewSet(current_pairIDs, allowed_pairs, curr_dimer_totals, nonset_dimers, None,
                                                  dimer_primerIDs, dimer_table, dimer_pairID, primer_pairs, primer_loci,
                                                  random=True, whitelist=whitelist_loci)

        # compare newSet to original set
        comparison, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total = compareSets(new_pairIDs, curr_total, swap_id, new_id, dimer_table, dimer_pairID)

        # if newSet is worse, make note of change value
        if comparison > 0:
            change.append(comparison)
            # repeat
            i += 1

    # B) set initial and final temps and temp schedule
    print("Setting temperature schedule...")
    # set initial and final temps
    T_end = min(change)
    # initial temp should accept most changes
    T_init = min(change) + ADJUSTMENT * (max(change) - min(change))
    T_init = int(round(T_init, 0))
    # partition iteration space across temp range
    if ITERATIONS > PARTITIONS:
        Tspace = (T_init - T_end)/PARTITIONS
        T_iter = Tspace*ITERATIONS
    else:
        Tspace = (T_init - T_end)/ITERATIONS
        T_iter = 1    
    print(".....Initial temp: "+str(T_init))
    print(".....Final temp: "+str(T_end))

    # C) run simulated annealing optimization
    print("Starting adaptive simulated annealing.....")
    Temp = T_init
    costs = [["Iterations","SA_Temp", "TotalDimers"]]
    costs.append([1, Temp, curr_total]) #log for cost values
    n_iter = 0
    step = 0
    min_dimers = curr_total
    while n_iter < ITERATIONS:
        i = 1  # reset iterations
        while i <= T_iter:
            # make a random change
            swap_id, new_id, new_pairIDs = MakeNewSet(current_pairIDs, allowed_pairs, curr_dimer_totals, nonset_dimers, None,
                                                      dimer_primerIDs, dimer_table, dimer_pairID, primer_pairs, primer_loci,
                                                      random=True, whitelist=whitelist_loci)
            # if e^(-change/temp) < random number (0,1), accept the change
            comparison, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total = compareSets(new_pairIDs, curr_total, swap_id, new_id, dimer_table, dimer_pairID)
            SAvalue = math.exp(-comparison/Temp)
            if SAvalue > rand.uniform(0,1):
                update = updateSet(swap_id, new_id, new_pairIDs, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total, verbose=False)
                current_pairIDs, curr_total, curr_dimer_totals, primerset_dimers, nonset_dimers = update #parse output into components
                # log updated cost whenever a change is made
                costs.append([n_iter+i, Temp, curr_total])
                # keep track of set with minimum dimers
                # this allows the best set to be kept while the algorithm continues to explore more of the space
                if curr_total < min_dimers:
                    min_pairIDs = current_pairIDs
                    min_dimers = curr_total
                    min_dimer_totals = curr_dimer_totals
                    min_primerset_dimers = primerset_dimers 
                    min_nonset_dimers = nonset_dimers
                # STOP if 0 dimers in set
                if curr_total == 0:
                    break
            # in simulated annealing, no need to keep track of what's been tested
            # repeat until # iterations at this temp finishes
            i+=1
        # decrease the temp via exponential decay - this makes the algorithm spend less time in 
        # 'risky' space and more time in productive optimization space
        step+=1
        Temp = (T_init-T_end)*DECAY_RATE**step+T_end
        n_iter+=T_iter
        #print("....."+str(int(n_iter))+" iterations")
        if curr_total == 0:
            print("Solution with 0 dimers found!")
            break
        print(curr_total)
    
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
    
    
    ### STEP 3: Simple iterative improvement with swapping worst values
    # only attempt iterative improvement if no final solution found with simulated annealing
    if curr_total > 0:
        print("Proceeding to simple iterative improvement....")
        print("Initial # dimers: " + str(curr_total))
        i = 0
        if len(whitelist_pairs) > 0:
            blacklist = whitelist_pairs
        else:
            blacklist = []
        while i < SIMPLE:
            # if current total dimers = 0, STOP because we're optimized already
            if curr_total == 0:
                print("STOPPING: Primer set is fully optimized.")
                break
    
            # create new set by replacing current worst pair with a random pair
            newSet = MakeNewSet(current_pairIDs, allowed_pairs, curr_dimer_totals, nonset_dimers, blacklist,
                                dimer_primerIDs, dimer_table, dimer_pairID, primer_pairs, primer_loci)
            if newSet is None:
                comparison = False
            else:
                swap_id, new_best_id, new_pairIDs = newSet
                # test if new set is better than current set (if fewer overall dimers)
                comparison, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total = compareSets(new_pairIDs, curr_total, swap_id, new_best_id, dimer_table, dimer_pairID)
    
            # if there are fewer dimers in new set than previous, then update the set
            if comparison < 0:
                update = updateSet(swap_id, new_id, new_pairIDs, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total, verbose=False)
                current_pairIDs, curr_total, curr_dimer_totals, primerset_dimers, nonset_dimers = update #parse output into components
                costs.append(['NA', curr_total])
    
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
                AddSwapLocusToAllowed(swap_id, allowed_pairs_rmv, primer_loci, primer_pairs)
                # avoid trying the same values...
                orig_swap = swap_id
                orig_best = new_best_id
                # loop through, allowing all other replacement options to be tested
                x = len(allowed_pairs_rmv)
                blacklist2 = []
                while x > 0:
                    if x > 0:
                        # set these to avoid infinite loops
                        prev_worst = swap_id
                        prev_best = new_best_id
                        # make new set
                        newSet = MakeNewSet(current_pairIDs, allowed_pairs_rmv, curr_dimer_totals, nonset_dimers, blacklist,
                                            dimer_primerIDs, dimer_table, dimer_pairID, primer_pairs, primer_loci, blacklist2)
                        if newSet is None:
                            print("No new sets can be found for "+swap_id+"! Blacklisting.")
                            blacklist.append(swap_id)
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
                        # add to blacklist2 if stuck in an infinite loop
                        if new_best_id == prev_best and swap_id == prev_worst:
                            print("Infinite loop caused by "+new_best_id+"! Blacklisting.")
                            blacklist2.append(new_best_id)
                        # compare new set against current set
                        comparison, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total = compareSets(new_pairIDs, curr_total, swap_id, new_best_id, dimer_table, dimer_pairID)
                        # keep new set if it has fewer dimers
                        if comparison < 0:
                            update = updateSet(swap_id, new_id, new_pairIDs, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, new_total, verbose=False)
                            current_pairIDs, curr_total, curr_dimer_totals, primerset_dimers, nonset_dimers = update #parse output into components
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
                    # when options run out for replacing this pair, add it to blacklist to avoid trying to replace
                    if x <= 1:
                        #print("Adding "+curr_worst+" to blacklist")
                        blacklist.append(swap_id)
    
            # stop iterating if no new sets can be made in the previous while loop
            if newSet is None:
                break
    
            # update iteration #
            i += 1
            if i==SIMPLE:
                print("Iterative improvement algorithm stopped because max # of iterations reached.")
            # progress tracking
            print(".....Current # dimers: "+str(curr_total))
    
        # report the blacklisted loci
        print("The following loci could not be improved with available replacements. If primer loads are high, consider removing these from the input file and rerunning.")
        for b in range(len(blacklist)):
            print("\t\t" + blacklist[b])
    
    
    ## STEP 4: Export selected primers, dimer loads, and cost trace
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
    with open(OUTPATH+'_primers.csv', 'w') as file:
        file.write("PrimerPairID,Sequence\n")
        for i in range(len(outpairs)):
            file.write(outpairs[i]+","+outseqs[i]+"\n")
    
    # export cost changes
    with open(OUTPATH+'_ASA_costs.csv', 'w') as file:
        for line in costs:
            file.write(str(line[0])+","+str(line[1])+","+str(line[2])+"\n")
    
    return curr_total



def MakeNewSet(pairIDs, allowed, curr_dimer_totals, nonset_dimers, blacklist, dimer_primerIDs, dimer_table,
               dimer_pairID, primer_pairs, primer_loci, blacklist2=[], random=False, whitelist=None):
    # create copy of original list, otherwise links them and affects both lists when using update/remove
    update_pairIDs = pairIDs.copy()
    allowed_pairs_edit = allowed.copy()

    # if swapping randomly, choose a random pairID (not including whitelist) to replace
    if random:
        options =  pairIDs.copy()
        for p in pairIDs:
            if p in whitelist:
                options.remove(p)
        swap = rand.choice(options)
    
    # otherwise, replace the current worst (most dimers) pair in the set
    else:
        worst = sorted(curr_dimer_totals.items(), key=lambda x: x[1])[-1]
        #curr_max = worst[1]
        swap = worst[0]

        # if the chosen ID is in the blacklist of loci that can't be fixed, then choose the next max value
        if swap in blacklist:
            j = -2
            n = len(curr_dimer_totals)
            while n >= abs(j):
                next_worst = sorted(curr_dimer_totals.items(), key=lambda x: x[1])[
                    j:(j+1)]  # get next worst
                #curr_max = next_worst[0][1]
                swap = next_worst[0][0]
                if swap not in blacklist:
                    break
                else:
                    j = j-1
            if n < abs(j):
                return None

    # add other primers for the swap locus back to the allowed list
    AddSwapLocusToAllowed(swap, allowed_pairs_edit, primer_loci, primer_pairs)

    # subset nonset dimer table to include only allowed primer pairs (i.e., primer pairs for loci not in the set)
    # gets all pair IDs that aren't in current pairs
    #nonset_pairIDs = list(set(dimer_primerIDs) - set(new_pairIDs))
    # gets indices for allowed pair IDs (i.e., not current loci) in nonset list
    #nonset_allowed_idx = list(filter(lambda x: dimer_primerIDs[x] in allowed_pairs_edit, range(len(dimer_primerIDs))))
    #nonset_pairs_allowed = [nonset_pairIDs[i] for i in nonset_allowed_idx]
    #nonset_dimers_sub = dict()
    # for pair in list(nonset_dimers.keys()):
    #    pairDict=nonset_dimers[pair]
    #    subPair=[pairDict[i] for i in nonset_allowed_idx]
    #    nonset_dimers_sub.update({pair: subPair})
    # ID primer pair outside of set that would have the smallest dimer load in the context of the current set
    # The main purpose of this is to speed up the optimization (as opposed to just choosing a random pair, many)
    # of which may not be very good in any context). This method also doesn't allow any flexibility if the min option doesn't work.
    # TODO: Fix this. Next set of lines gets stuck in an infinite loop without doing much optimization
    # nonset_dimers_sub.pop(curr_worst) #exclude pair being replaced from calculations
    # nonset_totals=list(map(sum, zip(*nonset_dimers_sub.values()))) #get sums across same column of values
    # nonset_min=min(nonset_totals)
    # grab the new pair ID
    #new_best_idx = list(filter(lambda x: nonset_totals[x]==nonset_min, range(len(nonset_totals))))
    # new_best_idx = random.choice(new_best_idx) #if more than 1 match, randomly choose 1...
    #new_best_id = nonset_pairs_allowed[new_best_idx]

    # Remove any primer pairs that cause infinite loops
    for pair in blacklist2:
        try:
            allowed_pairs_edit.remove(pair)
        except:
            pass
    # Choose random primer pair as replacement
    if len(allowed_pairs_edit) > 0:
        new_id = rand.choice(allowed_pairs_edit)
    else:
        return None
    # double check that chosen pair isn't already in set... if it is, remove and try again
    while new_id in update_pairIDs:
        allowed_pairs_edit.remove(new_id)
        if len(allowed_pairs_edit) > 0:
            new_id = rand.choice(allowed_pairs_edit)
        else:
            return None
    # update the primer set
    update_pairIDs.append(new_id)
    update_pairIDs.remove(swap)
    return swap, new_id, update_pairIDs


def compareSets(new_pairIDs, curr_total, swap, new_best_id, dimer_table, dimer_pairID):
    # calculate # dimers in this new set
    new_primerset_dimers, new_nonset_dimers = SubsetDimerTable(new_pairIDs, dimer_table, dimer_pairID, True)
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
    # else:
    #    return False


def AddSwapLocusToAllowed(swap_pair, inlist, primer_loci, primer_pairs):
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
    return inlist


def condition(x, list):
    if x in list:
        return x


def GetLocusID(pairID):
    pairIDsplit = pairID.split(".")
    LocusID = pairIDsplit[0]
    return LocusID


def BestPrimers(loci, dimer_loci, dimer_tallies, dimer_primerIDs, whitelist):
    """loci: """
    """dimer_tally_dict: """
    # replace "best" with whitelist pair if in whitelist
    whitelist_loci = [whitelist[x].split('.')[0] for x in range(len(whitelist))]
    best_primer_pairs = dict()
    for locus in loci:
        if locus in whitelist_loci:
            # assign whitelist pairs as best
            whitelist_indx = list(
                filter(lambda x: whitelist_loci[x] == locus, range(len(whitelist_loci))))
            bestPair = [whitelist[x] for x in whitelist_indx][0]
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


def SubsetDimerTable(primer_pairs, dimer_table, dimer_ids, return_complement=False):
    """primer_pairs"""
    """dimer_tally_dict"""
    sub_dimer_table = dict()  # blank dictionary to store output
    not_dimer_table = dict()
    # grab dictionary indices for the input primer pairs (skip first ID because that's the name column)
    primer_idx = list(
        filter(lambda x: dimer_ids[x] in primer_pairs, range(1, len(dimer_ids))))
    # left-adjust by 1 due to empty first column in dimer_ids
    primer_idx = [i-1 for i in primer_idx]
    for pair in primer_pairs:
        # grab all primer pair interactions between the specified pairs
        pair_dict = dimer_table[pair]
        sub_pair_dict = [int(pair_dict[i]) for i in primer_idx]
        sub_dimer_table.update({pair: sub_pair_dict})
        # grab all primer pair interactions between the pair and any pairs not specified
        not_pair_dict = [int(pair_dict[i])
                         for i in range(len(pair_dict)) if i not in primer_idx]
        not_dimer_table.update({pair: not_pair_dict})
    # return subset table (and table with non-subset values)
    if not return_complement:
        return sub_dimer_table
    else:
        return sub_dimer_table, not_dimer_table


def CalcTotalDimers(dimerDict):
    outDict = dict()
    for pair in dimerDict.keys():
        Ndimers = sum(dimerDict[pair])
        outDict.update({pair: Ndimers})
    return outDict


def LoadPrimers(PRIMER_FASTA, whitelist=False):
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
    if whitelist:
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


# run main function when called from the command line
if __name__ == "__main__":
    main(sys.argv[1],
         sys.argv[2],
         sys.argv[3],
         sys.argv[4],
         int(sys.argv[5]),
         int(sys.argv[6]),
         sys.argv[7])
