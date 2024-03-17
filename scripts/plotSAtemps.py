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
import matplotlib as mpl
import matplotlib.pyplot as plt
import random as rand
import csv
import sys
import math



def main(OUTPATH, PRIMER_FASTA=None, DIMER_SUMS=None, DIMER_TABLE=None, N_LOCI=None, WHITELIST=None, SEED=None, 
         MIN_DIMER=None, MAX_DIMER=None, DECAY_RATE=0.95, T_INIT=None, T_FINAL=None, BURNIN=100, ADJUSTMENT=0.1):
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
    ADJUSTMENT = Numeric (0-1)
        Proportion of max dimer load to consider when setting T_INIT (high values = more 'bad' changes accepted)
    -------
    Iteratively optimizes a multiplex primer set of given size to minimize
            predicted dimer content. NOTE: This script relies on having a large ratio of 
            loci to choose from relative to the number of loci needed.
    Outputs : CSVs of selected primer pairs + dimer loads
    """
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
        
            # read in whitelist info
            if WHITELIST is not None:
                whitelist_loci, whitelist_pairs = LoadPrimers(WHITELIST, True)
                n_whitelist = len(whitelist_pairs)
            else:
                whitelist_pairs = []
                whitelist_loci = []
                n_whitelist = 0
            
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
                                                          random=True, whitelist=whitelist_loci)
        
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
    
    ## Define range of dimer loads to use in calculations
    dimers = list(range(MIN_DIMER, MAX_DIMER+1))


    ## STEP 2: Set temperature schedule
    print("Setting temperature schedule...")
    if T_FINAL is None:
        # set initial and final temps
        T_FINAL = float(MIN_DIMER)
    if T_INIT is None:
        # initial temp should accept most changes
        T_INIT = MIN_DIMER + ADJUSTMENT * (MAX_DIMER - MIN_DIMER)
    # calculate temperatures across 100 iterations
    temps = [T_INIT]
    i=0
    while i < 100:
        i+=1
        new_temp = (T_INIT-T_FINAL)*DECAY_RATE**i+T_FINAL
        temps.append(new_temp)
    print(".....Initial temp for adaptive simulated annealing: "+str(T_INIT))
    print(".....Final temp for adaptive simulated annealing: "+str(T_FINAL))
    print("Plotting temperature schedule over 100 iterations")
    plt.plot(temps)
    plt.ylabel("Temperature")
    plt.xlabel("Iterations")
    plt.title("Temperature Schedule")
    plt.savefig(OUTPATH+"_TemperatureSchedule.png")
    
    
    ## STEP 4: Calculate and plot probabilities across range of dimer load changes
    print("Plotting acceptance probabilities for dimers across temperatures...")
    # grab color gradient
    reds = mpl.colormaps['Reds']
    cspace = 1 / len(dimers) # rescale dimer length to 0-1
    # plot each dimer load
    i = 0
    n = 0 # starting value for color scale
    while i<len(dimers):
        d = dimers[i]
        n+=cspace#update color scale value
        c = reds(n)
        plt.plot(temps, ([math.exp(-d/x) for x in temps]), color=c, label=str(d)+' dimers')
        i+=1
    # rescale axes and add labels
    plt.title("Dimer Acceptance as Simulated Annealing Progresses")
    plt.xlabel("Temperature decay")
    plt.ylabel("Probability of accepting dimers")
    plt.legend(loc='upper right')
    plt.xlim(T_INIT, T_FINAL)
    plt.ylim(0, 1)
    plt.savefig(OUTPATH+"_DimerAcceptanceProbs.png")



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



def GetLocusID(pairID):
    pairIDsplit = pairID.split(".")
    LocusID = pairIDsplit[0]
    return LocusID



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



if __name__=="__main__":
    main(sys.arv[1], sys.arv[2], sys.arv[3], sys.argv[4], sys.argv[5], sys.arv[6], sys.argv[7],
         sys.arv[8], sys.arv[9], sys.arv[10], sys.argv[11], sys.argv[12], sys.argv[13])