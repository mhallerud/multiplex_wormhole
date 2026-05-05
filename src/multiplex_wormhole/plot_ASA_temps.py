#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: PLOT SIMULATED ANNEALING PARAMETERS
Purpose: Assessing effects of simulated annealing parameter choices on acceptance
probability of dimer loads.
Dependencies: matplotlib (developed w/ v3.5.2)
              pandas (developed w/ v1.4.4)

Created on Sun Mar 17 13:21:40 2024
@author: maggiehallerud
"""


# load dependencies
import os
import sys
#import csv
import logging
from datetime import datetime
import math
import importlib
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import random as rand
import argparse


# load functions from optimize_primers module
sys.path.append(os.path.dirname(__file__))
op = importlib.import_module("optimize_multiplex")
from helpers.logging_setup import setup_logging



def main(OUTPATH, PRIMER_FASTA=None, DIMER_SUMS=None, DIMER_TABLE=None, N_LOCI=None, KEEPLIST=None, deltaG=False, SEED=None, 
         MIN_DIMER=None, MAX_DIMER=None, DECAY_RATE=0.95, T_INIT=None, T_FINAL=0.1, BURNIN=100, PROB_ADJ=2):
    """
    PRIMER_FASTA : Primer sequences and IDs. [FASTA]
    DIMER_SUMS : Total dimer load per primer pair. [CSV]
    DIMER_TABLE : Pairwise interaction table of primer pairs. [CSV]
    OUTPATH : Prefix for output files. [Filepath]
    N_LOCI : Number of primer pairs in final multiplex. [int]
    KEEPLIST : FASTA of primers required to be included in final multiplex [Default: None]
    deltaG : Optimize for total dimer tally (False) or maximum average deltaG among dimers (True)? [Default: False]
    SEED : Initial primer set to use as starting point, e.g., from previous run. [CSV] Default: None
    MIN_DIMER : Minimum # dimers for plotting acceptance probabilities. [int] Default: None
    MAX_DIMER : Maximum # dimers for plotting acceptance probabilities. [int] Default: None
    DECAY_RATE : Base for exponential decay of temperatures -- closer to 1 is slower decay rate. [num: 0-1; Default: 0.95]
    T_INIT : Starting simulated annealing temperature. [num>=0; Default: None]
    T_FINAL : Ending simulated annealing temperature. [num<T_INIT; Default: 0.1]
    BURNIN : Iterations used to sample cost space before calculating temperatures. [Default: 100]
    PROB_ADJ : Adjustment for dimer acceptance probabilities in exponential function (higher = fewer 'bad' changes accepted). [Default: 2]
    -------
    Plots temperature schedule and dimer acceptance probabilities based on given simulated annealing parameters.
    """
    ## first check files exist....
    CheckInput(PRIMER_FASTA, "PRIMER_FASTA")
    CheckInput(DIMER_SUMS, "DIMER_SUMS")
    CheckInput(DIMER_TABLE, "DIMER_TABLE")
    CheckInput(KEEPLIST, "KEEPLIST")
    CheckInput(SEED, "SEED")    
    
    # initialize logging
    logger = setup_logging(OUTPATH+".log", True, NAME=OUTPATH)    
    #print("Logging printed to: ")
    #for h in logger.handlers:
    #    print(type(h), getattr(h, OUTPATH+'.log', None))
    # log start time & inputs
    logger.info("START TIME: %s", datetime.now().strftime('%m/%d/%Y %I:%M:%S %p'))
    logger.info("")
    logger.info("plot_ASA_temps inputs: ")
    logger.info("     OUTPATH: %s", OUTPATH)
    logger.info("     PRIMER_FASTA: %s", PRIMER_FASTA)
    logger.info("     DIMER_SUMS: %s", DIMER_SUMS)
    logger.info("     DIMER_TABLE: %s", DIMER_TABLE)
    logger.info("     N_LOCI: %s", N_LOCI)
    logger.info("     KEEPLIST: %s", KEEPLIST)
    logger.info("     deltaG: %s", deltaG)
    logger.info("     SEED: %s", SEED)
    logger.info("     MIN_DIMER: %s", MIN_DIMER)
    logger.info("     MAX_DIMER: %s", MAX_DIMER)
    logger.info("     DECAY_RATE: %s", DECAY_RATE)
    logger.info("     BURNIN: %s", BURNIN)
    logger.info("     T_INIT: %s", T_INIT)
    logger.info("     T_FINAL: %s", T_FINAL)
    logger.info("     PROB_ADJ: %s", PROB_ADJ)
    logger.info("")
    
    ## OPTION a: SIMULATED ANNEALING WITH FIXED SCHEDULE
    if T_INIT is not None and T_FINAL is not None:
        pass
    
    ## OPTION B: SIMULATED ANNEALING WITH ADAPTIVE SCHEDULE
    else:
        ## STEP 1: Calculate min and max observed dimer load changes (if not provided)
        if MAX_DIMER is None and MIN_DIMER is None:
            # check for necessary files before proceeding...
            if not [PRIMER_FASTA is not None and DIMER_SUMS is not None and DIMER_TABLE is not None and N_LOCI is not None]:
                raise InputError("If MAX_DIMER and MIN_DIMER are not specified, paths to PRIMER_FASTA, DIMER_SUMS, "+
                                "DIMER_TABLE and N_LOCI must be provided. Please provide either A) MAX_DIMER and MIN_DIMER or B)"+
                                " PRIMER_FASTA, DIMER_SUMS, DIMER_TABLE, and N_LOCI parameters.")
            else:
                ## Read in IDs and primers
                logger.info("Reading in inputs..........")
                # read in IDs and primers
                primer_loci, primer_seqs, primer_IDs, primer_pairs = op.LoadPrimers(PRIMER_FASTA)

            
                # read in dimer info
                #logger.info("Reading in primer dimer counts........")
                dimer_table = pd.read_csv(DIMER_TABLE)
                dimer_sums = pd.read_csv(DIMER_SUMS)
            
                # read in keeplist info
                if KEEPLIST is not None:
                    keeplist_loci, keeplist_seqs, keeplist_IDs, keeplist_pairs = op.LoadPrimers(KEEPLIST)
                    n_keeplist = len(set(keeplist_pairs))
                else:
                    keeplist_pairs = []
                    #keeplist_loci = []
                    keeplist_seqs = []
                    keeplist_IDs = []
                    n_keeplist = 0
                
                # check for missing dimer info
                missing_pairs = [keeplist_pairs[x] for x in range(len(keeplist_pairs)) if keeplist_pairs[x] not in dimer_table['Pair1'].to_numpy()]
                if len(missing_pairs)>0:
                    raise InputError("Dimer information is missing for keeplist pairs! Rerun MFEprimer dimer prediction after running AddKeeplist2FASTA")
                primer_pairs, primer_loci, primer_seqs, primer_IDs = op.CheckPrimersInDimerTables(
                    dimer_table, primer_pairs, primer_loci, primer_seqs, primer_IDs, logger)

                if SEED is None:
                    ## A) Choose initial primer set based on pseudo-Greedy algorithm
                    # Choose initial set of loci based on loci with minimum dimer counts
                    # (Hopefully choosing an initial set this way, rather than randomly, will mean fewer iterations are needed)
                    logger.info("Generating initial primer set........")
                    # make list of unique loci
                    # convert to list because new versions of random.sample won't be able to handle sets...
                    uniq_loci = list(set(primer_loci))
                    nloci = len(uniq_loci)
                    # grab best primer pairs for each locus
                    best_primer_pairs = op.BestPrimers(uniq_loci, dimer_sums, keeplist_pairs, deltaG)
                    # if there are fewer loci than desired, use all of them
                    if nloci < N_LOCI:
                        logger.info("WARNING: Fewer loci passed filtering than desired in panel")
                        logger.info("# loci used: %s", str(nloci))
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
                    logger.info("Loading initial primer set from seed file........")
                    initial_pairs = []
                    with open(SEED, 'r') as file:
                        lines = file.readlines()
                        for line in lines:
                            initial_pairs.append(lines[0])
            
                # grab locus IDs for primer pairs
                current_pairIDs = list(initial_pairs.keys())
                current_pairIDs = [str(x) for x in current_pairIDs]
                current_locusIDs = [op.GetLocusID(pair) for pair in current_pairIDs]
                # get initial list of allowed alternative primer pairs (i.e., primer pairs for loci not currently in set)
                allowed_loci = [uniq_loci[i] for i in range(
                    len(uniq_loci)) if uniq_loci[i] not in current_locusIDs]
                allowed_idx = list(
                    filter(lambda x: primer_loci[x] in allowed_loci, range(len(primer_loci))))
                allowed_pairs = [primer_pairs[i] for i in allowed_idx]
                allowed_pairs = list(set(allowed_pairs))
                # calculate # primer dimers per pair for current set of primers
                primerset_dimers, nonset_dimers = op.SubsetDimerTable(current_pairIDs, dimer_table, True)
                curr_dimer_totals = op.CalcTotalDimers(primerset_dimers, deltaG)  
                curr_total = sum(curr_dimer_totals.values())
                
                
                ## B) Sample cost space
                if nloci <= N_LOCI:
                    logger.info("FEWER LOCI INPUT THAT DESIRED IN FINAL PANEL: "+
                          "Proceeding to plot with default temps with 1-5 dimers, "+
                          "but true optimization for this problem set would stop prior to the simulated annealing step.")
                    T_INIT = 2
                    T_FINAL = 0
                    MIN_DIMER = 1
                    MAX_DIMER = 5
                else:
                    logger.info("Sampling dimer load changes.......")
                    # A) sample x iterations
                    change = []
                    i = 0
                    while len(change) < BURNIN:
                        # make a new set by randomly swapping a primer pair
                        # newset: 1) replaced ID, 2) new ID, 3) current pair list
                        swap_id, new_id, new_pairIDs = MakeNewSet(current_pairIDs, allowed_pairs, curr_dimer_totals, nonset_dimers, [],
                                                                  primer_pairs, primer_loci, 
                                                                  OUTPATH, primer_IDs, primer_seqs, keeplist_IDs, keeplist_seqs, [], [],
                                                                  random=True, keeplist=keeplist_pairs)
                        # compare newSet to original set
                        comparison, new_primerset_dimers, new_nonset_dimers, new_dimer_totals, \
                            new_total = op.compareSets(new_pairIDs, curr_total, swap_id, 
                                                    new_id, dimer_table, dimer_sums, deltaG)
                
                        # if newSet is worse, make note of change value
                        if comparison > 0:
                            change.append(comparison)
                            
                        # stop if 2000 swaps made
                        i += 1
                        if i>2000:
                            logger.info("....Sampling stopped after 2000 swaps with %s costs>0.", str(len(change)))
                            if len(change)==0:
                                logger.info(".....Defaults used since cost changes never increase")
                                T_INIT = 2
                                T_FINAL = 0
                                MIN_DIMER=1
                                MAX_DIMER=5
                            break
                    if len(change)>0:
                        MIN_DIMER = min(change)
                        MAX_DIMER = max(change)
                        logger.info(".....Maximum dimer load observed %s", str(MAX_DIMER))
                        logger.info(".....Minimum dimer load observed %s", str(MIN_DIMER))


        ## STEP 2: Set adaptive temperature schedule based on cost calcs
        logger.info("Setting temperature schedule...")
        if T_FINAL is None:
            # set initial and final temps
            T_FINAL = 0.01
        if T_INIT is None:
           T_INIT, MIN_DIMER, MAX_DIMERS =  op.setTemps(current_pairIDs, allowed_pairs, curr_dimer_totals, nonset_dimers, 
                                 primer_pairs, primer_loci, OUTPATH, primer_IDs, primer_seqs, keeplist_IDs, 
                                 keeplist_seqs, keeplist_pairs=keeplist_pairs, RNG=12345, 
                                 curr_total=curr_total, dimer_table=dimer_table, dimer_sums=dimer_sums, 
                                 deltaG=deltaG, logger=logger)
        else:
             T_INIT = float(T_INIT)
    
    ## PLOT 1: GENERATE TEMPERATURE SCHEDULE ACROSS 100 steps
    # Define range of dimer loads to use in calculations
    dimers = list(range(round(MIN_DIMER), round(MAX_DIMER)+1))
    temps = [T_INIT]
    step=0
    while step < 100:
        step+=1
        new_temp = (T_INIT-T_FINAL)*DECAY_RATE**step+T_FINAL
        temps.append(new_temp)
    logger.info(".....Initial temp for adaptive simulated annealing: %s", str(T_INIT))
    logger.info(".....Final temp for adaptive simulated annealing: %s", str(T_FINAL))
    logger.info("Plotting temperature schedule over 100 iterations")
    plt.figure() # new plotting figure
    plt.plot(temps)
    plt.ylabel("Temperature")
    plt.xlabel("% of Iterations")
    plt.title("Temperature Schedule")
    plt.savefig(OUTPATH+"_TemperatureSchedule.png")
    
    
    ## PLOT 2: DIMER ACCEPTANCE PROBABILITY BY SA TEMPERATURE
    logger.info("Plotting acceptance probabilities for dimers across temperatures...")
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
        costs.append([math.exp(-PROB_ADJ*d/x) for x in temps])
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
    
    # end logging
    logger.info("")
    logger.info("END TIME: %s", datetime.now().strftime('%m/%d/%Y %I:%M:%S %p'))
    logging.shutdown()



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
    #worst_dimers = [x[1] for x in worst]
    #worst_sum = sum(worst_dimers)
    #if worst_sum==0:
    #    costs.append([n_iter, Temp, curr_total])
    #    ExportCSVs(OUTPATH, curr_dimer_totals, primer_pairs, pairIDs, primer_IDs, 
    #               primer_seqs, keeplist_IDs, keeplist_seqs, costs)
    #    raise OptimizationWarning("STOPPED. Keeplist primer pairs are the only ones with dimer loads, therefore further optimization is impossible. (outputs saved).")

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
    allowed_pairs_edit = op.UpdateAllowedPairs(swap, allowed_pairs_edit, primer_loci, primer_pairs)
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



def CheckInput(file, name, required=False):
    if required or file is not None:
        if not os.path.exists(file):
            raise InputError(name+" file could not be found!")



class InputError(Exception):
    pass



def parse_args():
    # initialize argparser
    parser = argparse.ArgumentParser()
    # add required arguments
    parser.add_argument("-o", "--outpath", type=str, required=True)
    # add optional arguments
    parser.add_argument("-f", "--primer_fasta", type=str, default=None)
    parser.add_argument("-d", "--dimer_sums", type=str, default=None)
    parser.add_argument("-t", "--dimer_table", type=str, default=None)
    parser.add_argument("-n", "--nloci", type=int, default=None)
    parser.add_argument("-k", "--keeplist", type=str, default=None)
    parser.add_argument("-e", "--seed", type=str, default=None)
    parser.add_argument("-m", "--min_dimers", type=int, default=None)
    parser.add_argument("-a", "--max_dimers", type=int, default=None)
    parser.add_argument("-r", "--decay_rate", type=float, default=0.95)
    parser.add_argument("-i", "--temp_init", type=float, default=None)
    parser.add_argument("-j", "--temp_final", type=float, default=0.1)
    parser.add_argument("-b", "--burnin", type=int, default=100)
    parser.add_argument("-p", "--prob_adj", type=float, default=2)
    # add flags
    parser.add_argument("-g", "--deltaG", action="store_true")
    return parser.parse_args()



if __name__=="__main__":
    # parse command line arguments
    args = parse_args()
    # run main
    main(OUTPATH = args.outpath, 
         PRIMER_FASTA = args.primer_fasta, 
         DIMER_SUMS = args.dimer_sums, 
         DIMER_TABLE = args.dimer_table, 
         N_LOCI = args.nloci, 
         KEEPLIST = args.keeplist,
         deltaG = args.deltaG,
         SEED = args.seed, 
         MIN_DIMER = args.min_dimers, 
         MAX_DIMER = args.max_dimers, 
         DECAY_RATE = args.decay_rate, 
         T_INIT = args.temp_init, 
         T_FINAL = args.temp_final, 
         BURNIN = args.burnin, 
         PROB_ADJ = args.prob_adj)
