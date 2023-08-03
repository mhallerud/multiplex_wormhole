#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 07:59:12 2023

@author: maggiehallerud
"""

# load dependencies
#import os
import sys
import random


# set user-defined parameters
PRIMER_FASTA = sys.argv[1]
DIMER_SUMS =sys.argv[2]
DIMER_TABLE = sys.argv[3]
N_LOCI = int(sys.argv[4])
ITERATIONS = int(sys.argv[5])

#PRIMER_FASTA="Test1_adapters" #filepath to adapter-ligated primers
#DIMER_SUMS="PrimerDimerReport_25jul2023_PrimerPairInteractions_sum_binary.csv"
#DIMER_TABLE="PrimerDimerReport_25jul2023_PrimerPairInteractions_wide_binary.csv"
#N_LOCI=10 # Number of SNPs in final panel
#ITERATIONS=50 # Number of times to run optimization



def main():
    ## read in IDs and primers
    ## NOTE: Numpy arrays would be more efficient & cleaner coding here than lists, 
    ## however numpy is constantly being updated so using it as a dependency would 
    ## make the script liable to break in the future. The functions herein are not 
    ## doing extremely heavy lifting, so indexing lists should work and provide the 
    ## most generalizable / user-friendly function.

    print("Reading in primer sequences..........")
    primer_loci, primer_seqs, primer_IDs = LoadPrimers(PRIMER_FASTA)
    
    ## read in dimer info
    print("Reading in primer dimer counts........")
    dimer_table, dimer_primerIDs, dimer_loci, dimer_tallies, dimer_pairID = LoadDimers(DIMER_SUMS, DIMER_TABLE)    
    
    ## Choose initial set of loci based on loci with minimum dimer counts
    ## (Hopefully choosing an initial set this way, rather than randomly, will mean fewer iterations are needed)
    print("Generating initial primer set........")
    #TODO: lock in important ones (e.g. microhaplotpyes, sex ID)
    #TODO: add ability to count SNPs within microhaps separately...
    # make list of unique loci
    uniq_loci=list(set(primer_loci)) # convert to list because new versions of random.sample won't be able to handle sets... 
    nloci=len(uniq_loci)
    # grab best primer pairs for each locus
    best_primer_pairs = BestPrimers(uniq_loci, dimer_loci, dimer_tallies, dimer_primerIDs)
    # if there are fewer loci than desired, use all of them
    if nloci < N_LOCI:
        print("WARNING: Fewer loci passed filtering than desired in panel")
        print("# loci used: " + str(nloci))
        initial_pairs = best_primer_pairs
    else:
        # grab N primer pairs with min dimer count
        initial_pairs = dict(sorted(best_primer_pairs.items(), key = lambda x: x[1])[:N_LOCI])
        #initial_loci=random.sample(uniq_loci, N_LOCI)
    
    # grab locus IDs for primer pairs
    current_pairIDs = list(initial_pairs.keys())
    current_locusIDs = [GetLocusID(pair) for pair in current_pairIDs]   
    # get initial list of allowed alternative primer pairs (i.e., primer pairs for loci not currently in set)
    allowed_loci = [uniq_loci[i] for i in range(len(uniq_loci)) if uniq_loci[i] not in current_locusIDs]
    allowed_idx = list(filter(lambda x: primer_loci[x] in allowed_loci, range(len(primer_loci))))
    allowed_pairs = [primer_IDs[i] for i in allowed_idx]
    allowed_pairs = list(set(allowed_pairs))
    # calculate # primer dimers per pair for current set of primers
    primerset_dimers, nonset_dimers = SubsetDimerTable(current_pairIDs, dimer_table, dimer_pairID, True)
    curr_dimer_totals = CalcTotalDimers(primerset_dimers)# totals for current primers
    curr_total = sum(curr_dimer_totals.values())
    
    ## Optimize primer set iteratively
    i=0
    blacklist=[]
    print("Starting optimization process......")
    print("Current # dimers: " + str(curr_total))
    while i < ITERATIONS:
        # if current total dimers = 0, STOP because we're optimized already
        if curr_total==0:
            print("STOPPING: Primer set is fully optimized.")
            #break

        # create new set by replacing current worst pair with a random pair
        newSet = MakeNewSet(current_pairIDs, allowed_pairs, curr_dimer_totals, nonset_dimers, blacklist, dimer_primerIDs, dimer_table, dimer_pairID, primer_IDs, primer_loci)
        if newSet is None:
            comparison=False
        else:
            curr_worst, new_best_id, new_pairIDs = newSet
            # test if new set is better than current set (if fewer overall dimers)
            comparison = compareSets(new_pairIDs, curr_total, curr_worst, new_best_id, allowed_pairs, dimer_table, dimer_pairID)

        #if test passed, then new set stats were returned
        if comparison is not False:
            current_pairIDs, curr_total, curr_dimer_totals, primerset_dimers, nonset_dimers = comparison
        
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
            AddWorstToAllowed(curr_worst, allowed_pairs_rmv, primer_loci, primer_IDs)
            # avoid trying the same values...
            orig_worst = curr_worst
            orig_best = new_best_id
            # loop through, allowing all other replacement options to be tested
            x=len(allowed_pairs_rmv)
            while x>0:
                if x>0:
                    # make new set
                    newSet = MakeNewSet(current_pairIDs, allowed_pairs_rmv, curr_dimer_totals, nonset_dimers, blacklist, dimer_primerIDs, dimer_table, dimer_pairID, primer_IDs, primer_loci)
                    if newSet is None:
                        print("No new sets can be found for this set of primer pairs")
                        break
                    else:
                        curr_worst, new_best_id, new_pairIDs = newSet
                    # skip back to top of while if new best is same as original curr_worst
                    if new_best_id==orig_worst or new_best_id==orig_best:
                        if x>1:
                            continue#go back to while
                        else:
                            break
                    # compare against current set
                    comparison = compareSets(new_pairIDs, curr_total, curr_worst, new_best_id, allowed_pairs, dimer_table, dimer_pairID)
                    if comparison is not False:
                        current_pairIDs, curr_total, curr_dimer_totals, primerset_dimers, nonset_dimers = comparison
                        break #exit this loop if a replacement was found
                    else:
                        # remove tested pair from allowed list
                        try:
                            allowed_pairs_rmv.remove(new_best_id)
                        except:
                            pass
                # otherwise, loop through remaining replacement options
                # this avoids retesting the same options for this replacement
                x=len(allowed_pairs_rmv)#update
                # when options run out for replacing this pair, add it to blacklist to avoid trying to replace
                if x==0:
                    #print("Adding "+curr_worst+" to blacklist")
                    blacklist.append(curr_worst)
        
        # stop iterating if no new sets can be made in the previous while loop
        if newSet is None:
            break

        # update iteration #
        i+=1
        # progress tracking
        print("\tCurrent # dimers: "+str(curr_total))

    # report the blacklisted loci
    print("The following loci could not be improved with existing replacements. If primer loads are high, consider removing these from the list.")
    for b in range(len(blacklist)):
        print("\t\t" + blacklist[b])
    
    return curr_dimer_totals



def MakeNewSet(pairIDs, allowed_pairs, curr_dimer_totals, nonset_dimers, blacklist,
               dimer_primerIDs, dimer_table, dimer_pairID, primer_IDs, primer_loci):
        # create copy of original list, otherwise links them and affects both lists when using update/remove    
        new_pairIDs = pairIDs.copy()
        allowed_pairs_edit = allowed_pairs.copy()
        # subset nonset dimer table to include only allowed primer pairs
        nonset_pairIDs = list(set(dimer_primerIDs) - set(new_pairIDs)) #gets pair IDs that aren't in current pairs
        nonset_allowed_idx = list(filter(lambda x: nonset_pairIDs[x] in allowed_pairs_edit, range(len(nonset_pairIDs))))
        nonset_pairs_allowed = [nonset_pairIDs[i] for i in nonset_allowed_idx]
        nonset_dimers_sub = dict()
        for pair in list(nonset_dimers.keys()):
            pairDict=nonset_dimers[pair]
            subPair=[pairDict[i] for i in nonset_allowed_idx]
            nonset_dimers_sub.update({pair: subPair})
        # ID current worst primer pair in set
        worst = sorted(curr_dimer_totals.items(), key = lambda x: x[1])[-1]
        curr_max = worst[1]
        curr_worst = worst[0]
        # if the chosen ID is in the blacklist of loci that can't be fixed, then choose the next max value
        if curr_worst in blacklist:
            j=-2
            n=len(curr_dimer_totals)
            while n>=abs(j):
                next_worst = sorted(curr_dimer_totals.items(), key = lambda x: x[1])[j:(j+1)]#get next worst
                curr_max = next_worst[0][1]
                curr_worst = next_worst[0][0]
                if curr_worst not in blacklist:
                    break
                else:
                    j=j-1
            if n<abs(j):
                return None
        # add primer IDs for the chosen locus to the allowed list
        AddWorstToAllowed(curr_worst, allowed_pairs_edit, primer_loci, primer_IDs)
        # ID primer pair outside of set that would have the smallest dimer load in the context of the current set
        # The main purpose of this is to speed up the optimization (as opposed to just choosing a random pair, many)
        # of which may not be very good in any context). This method also doesn't allow any flexibility if the min option doesn't work.
        #TODO: Fix this. Next set of lines gets stuck in an infinite loop without doing much optimization
        #nonset_dimers_sub.pop(curr_worst) #exclude pair being replaced from calculations
        #nonset_totals=list(map(sum, zip(*nonset_dimers_sub.values()))) #get sums across same column of values
        #nonset_min=min(nonset_totals)
        # grab the new pair ID
        #new_best_idx = list(filter(lambda x: nonset_totals[x]==nonset_min, range(len(nonset_totals))))
        #new_best_idx = random.choice(new_best_idx) #if more than 1 match, randomly choose 1...
        #new_best_id = nonset_pairs_allowed[new_best_idx]
        # Choose random primer pair as replacement
        new_best_id = random.choice(allowed_pairs_edit)
        # update the primer set
        new_pairIDs.append(new_best_id)
        new_pairIDs.remove(curr_worst)
        return curr_worst, new_best_id, new_pairIDs


def compareSets(new_pairIDs, curr_total, curr_worst, new_best_id, allowed_pairs, dimer_table, dimer_pairID):
    # calculate # dimers in this new set
    new_primerset_dimers, new_nonset_dimers = SubsetDimerTable(new_pairIDs, dimer_table, dimer_pairID, True)
    new_dimer_totals = CalcTotalDimers(new_primerset_dimers)
    new_total = sum(new_dimer_totals.values())
    # if this primer set has fewer dimers overall, keep it and proceed
    if new_total < curr_total:
        print("\t\t"+curr_worst + " replaced with " + new_best_id)
        # TODO: if a primer pair has been tried, remove from current 'allowed' set (so that we don't go in circles...)
        # update all starting values to match this iteration's values
        current_pairIDs = new_pairIDs.copy()
        try:
            allowed_pairs.remove(new_best_id)
        except:
            pass
        #allowed_pairs.append(curr_worst)
        primerset_dimers, nonset_dimers = new_primerset_dimers, new_nonset_dimers
        curr_dimer_totals = new_dimer_totals
        curr_total = new_total
        return current_pairIDs, curr_total, curr_dimer_totals, primerset_dimers, nonset_dimers
    else:
        return False


def AddWorstToAllowed(curr_worst, inlist, primer_loci, primer_IDs):
        curr_worst_locus = curr_worst.split("_")[0] + "_" + curr_worst.split("_")[1]
        curr_worst_idx = list(filter(lambda x: primer_loci[x]==curr_worst_locus, range(len(primer_loci))))
        curr_worst_pairs = [primer_IDs[i] for i in curr_worst_idx]
        curr_worst_pairs = (set(curr_worst_pairs))#extract unique values only
        curr_worst_pairs.remove(curr_worst)
        for p in curr_worst_pairs:
            inlist.append(p)


def condition(x, list): 
    if x in list:
        return x


def GetLocusID(pairID):
    pairIDsplit = pairID.split("_")
    LocusID = pairIDsplit[0] + "_" + pairIDsplit[1]
    return LocusID


def BestPrimers(loci, dimer_loci, dimer_tallies, dimer_primerIDs):
    """loci: """
    """dimer_tally_dict: """
    best_primer_pairs= dict()
    for locus in loci:
        # extract dimer data for this locus
        locus_idx = list(filter(lambda x: dimer_loci[x]==locus, range(len(dimer_loci))))
        ids = [dimer_primerIDs[i] for i in locus_idx]
        tallies = [dimer_tallies[i] for i in locus_idx]
        # extract primer pair with min dimer value
        min_idx = list(filter(lambda x: tallies[x]==min(tallies), range(len(tallies))))
        if len(min_idx)>1:# if there's more than one 'best' option, randomly choose one...
            min_idx = [random.choice(min_idx)]
        # extract the ID for the selected pair
        bestPair= [ids[i] for i in min_idx]
        # set up dictionary with locus: (primer pair id, min # dimers)
        best_primer_pairs.update({bestPair[0]: min(tallies)})
    return best_primer_pairs


def SubsetDimerTable(primer_pairs, dimer_table, dimer_ids, return_complement=False):
    """primer_pairs"""
    """dimer_tally_dict"""
    sub_dimer_table = dict() #blank dictionary to store output
    not_dimer_table = dict()
    # grab dictionary indices for the input primer pairs (skip first ID because that's the name column)
    primer_idx = list(filter(lambda x: dimer_ids[x] in primer_pairs, range(1,len(dimer_ids))))
    primer_idx = [i-1 for i in primer_idx]#left-adjust by 1 due to empty first column in dimer_ids
    for pair in primer_pairs:
        # grab all primer pair interactions between the specified pairs
        pair_dict = dimer_table[pair]
        sub_pair_dict = [int(pair_dict[i]) for i in primer_idx]
        sub_dimer_table.update({pair: sub_pair_dict})
        # grab all primer pair interactions between the pair and any pairs not specified
        not_pair_dict = [int(pair_dict[i]) for i in range(len(pair_dict)) if i not in primer_idx]
        not_dimer_table.update({pair: not_pair_dict})
    # return subset table (and table with non-subset values)
    if not return_complement:
        return sub_dimer_table
    else:
        return sub_dimer_table, not_dimer_table
    

def CalcTotalDimers(dimerDict):
    outDict=dict()
    for pair in dimerDict.keys():
        Ndimers=sum(dimerDict[pair])
        outDict.update({pair: Ndimers})
    return outDict


def LoadPrimers(PRIMER_FASTA):
    # read in files as lists
    primer_loci = []
    primer_IDs = []
    primer_seqs=[]
    # parse sequences from headers
    with open(PRIMER_FASTA, "r") as file:
        lines = file.readlines()
        for line in lines:
            if ">" in line:
                line=line.rstrip()
                line=line.replace('>','')
                linesplit=line.split("_")
                locus=linesplit[0]+"_"+linesplit[1]
                primer_loci.append(locus)
                primer_IDs.append(locus+"_"+linesplit[2])
            else:
                line=line.rstrip()
                primer_seqs.append(line)
    return primer_loci, primer_seqs, primer_IDs


def LoadDimers(DIMER_SUMS, DIMER_TABLE):
    dimer_primerIDs=[]
    dimer_loci=[]
    dimer_tallies=[]
    with open(DIMER_SUMS, 'r') as file:
        lines = file.readlines()
        for line in lines[1:]: #skip first line because R puts unnecessary header
            # read in line and reformat
            line=line.rstrip()
            linesplit=line.split(",")
            linesplit[0]=linesplit[0].replace('"', '')#remove extra quotes from string
            linesplitid=linesplit[0].split("_")
            # extract info
            dimer_primerIDs.append(linesplit[0])
            dimer_loci.append(linesplitid[0]+"_"+linesplitid[1])
            dimer_tallies.append(int(linesplit[1]))
    # we'll store this in a dictionary so that it's really easy to subset based on the primer pair ID
    dimer_table = dict()
    dimer_pairID = []
    with open(DIMER_TABLE, 'r') as file:
        lines = file.readlines()
        for line in lines:
            line=line.rstrip()
            linesplit=line.split(",")
            linesplit[0]=linesplit[0].replace('"', '')
            dimer_table.update({linesplit[0]: linesplit[1:]})
            dimer_pairID.append(linesplit[0])
    return dimer_table, dimer_primerIDs, dimer_loci, dimer_tallies, dimer_pairID



# run main function when called from the command line
if __name__ == "__main__":
    main()