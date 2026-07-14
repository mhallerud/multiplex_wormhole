#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: MULTIPLE OPTIMIZATION RUN
Purpose: Runs optimize_multiplex multiple times, organizes outputs and produces a summary.
    
Created on Fri Mar 15 14:31:28 2024
@author: maggiehallerud
"""

import traceback
import os
import sys
import glob
import shutil
import argparse
from datetime import datetime
import csv
from multiprocessing.pool import Pool
from functools import partial
from math import ceil

try:
    from .optimize_multiplex import main as optimizeMultiplex
    from .optimize_multiplex import CheckInputFile
    from .optimize_multiplex import LoadPrimers
    from .panel_assessment import main as assessPanel
except ImportError:
    sys.path.append(os.path.dirname(__file__))
    from optimize_multiplex import main as optimizeMultiplex
    from optimize_multiplex import CheckInputFile
    from optimize_multiplex import LoadPrimers
    from panel_assessment import main as assessPanel



def main(N_RUNS, PRIMER_FA, DIMER_SUMS, DIMER_TABLE, OUTPATH, N_LOCI, 
         deltaG=False, KEEPLIST=None, VERBOSE=False, SEED=None, #TIMEOUT=5,
         SIMPLE=5000, ITERATIONS=1000, CYCLES=10, BURNIN=200, DECAY_RATE=0.95, 
         T_INIT=None, T_FINAL=0.001, PROB_ADJ=2, MAKEPLOT=False,
         THREADS=None, dG_END_LIMIT=-4, dG_MID_LIMIT=-8, dG_BAD_LIMIT=-10):
    """
    N_RUNS : # optimization runs [int]
    PRIMER_FA : Contains primer IDs and sequences [FASTA]
    DIMER_SUMS : Sum of interactions per primer pair (binary interactions recommended) [CSV]
    DIMER_TABLE : Pairwise interaction table of primer pairs (binary interactions recommended) [CSV]
    OUTPATH : Output path and prefix [filepath]
    N_LOCI : Number of primer pairs in desired multipled [integer]
    KEEPLIST : Contains primers that MUST be included in the final solution [FASTA; default=None]
    deltaG : Minimize mean deltaG [True] or count of dimers- requires deltaG dimer tables! [Default: False]
    VERBOSE : Print progress? [Default: False]
    SEED : Initial primer set to start with from previous multiplex_wormhole run [CSV]
    -------SIMULATED ANNEALING PARAMETERS-----
    SIMPLE : # iterations in simple iterative improvement optimization step [default=5000]
    ITERATIONS : # iterations per simulated annealing optimization cycle [default=1000]
    CYCLES : # simulated annealing cycles to run [default: 10]
    BURNIN : # iterations to sample dimer cost space [integer; default=200]
    DECAY_RATE : parameter for temperature decay rate in negative exponential [proportion; default=0.95]
        -closer to 1: least conservative, explores more of cost space but adds more dimers
    T_INIT : Starting temperature for simulated annealing (default: None, i.e. set adaptively)
        - values closer to 0 will reduce increasing costs
    T_FINAL : Ending temperature for simulated annealing (default: 0.001)
    PROB_ADJ : Parameter used to adjust dimer acceptance probabilities (default: 2)
        -Try increasing to 2 or 3 if too many dimers are being accepted during simulated annealing
    -------
    Outputs N_RUNS optimized multiplexes and a summary.

    """
    # first check that files exist
    CheckInputFile(PRIMER_FA, "PRIMER_FA", required=True)
    CheckInputFile(DIMER_SUMS, "DIMER_SUMS", required=True)
    CheckInputFile(DIMER_TABLE, "DIMER_TABLE", required=True)
    CheckInputFile(KEEPLIST, "KEEPLIST", required=False)
    CheckInputFile(SEED, "SEED", required=False)
    
    # check that primer IDs are in correct format
    try:
        test = LoadPrimers(PRIMER_FA)
    except Exception:
        raise InputError("PrimerIDs in PRIMER_FA are not in the correct format! " \
                         "Reformat IDs as <name>.<#>.<DIR>, e.g., MACA01.0.FWD & MACA01.0.REV")
    
    if KEEPLIST is not None:
        try:
            test = LoadPrimers(KEEPLIST)
        except Exception:
            raise InputError("PrimerIDs in KEEPLIST are not in the correct format! " \
                             "Reformat IDs as <name>.<#>.<DIR>, e.g., MACA01.0.FWD & MACA01.0.REV")
    
    if THREADS>N_RUNS:
        print("More THREADS specified than N_RUNS - CPU usage will be limited to the # of runs.")
    # set up output filename if None
    if OUTPATH is None or OUTPATH=="None":
        OUTPATH = datetime.now().strftime("%d-%m-%Y_%H%M")
    # set up empty array to hold overall dimer load 
    loads = [['Run', 'Filepath','N_Pairs','TotalDimers','PairsWithDimers','DimersPerPair',
              'MeanDeltaG','N_BadDimers']]
        
    # run optimizations with multi-processing
    # set up arguments as list
    worker = partial(runOpt, PRIMER_FA=PRIMER_FA, DIMER_SUMS=DIMER_SUMS, DIMER_TABLE=DIMER_TABLE,
                     OUTPATH=OUTPATH, N_LOCI=N_LOCI, KEEPLIST=KEEPLIST, deltaG=deltaG, SEED=SEED,
                     VERBOSE=VERBOSE, SIMPLE=SIMPLE, ITERATIONS=ITERATIONS, CYCLES=CYCLES,
                     BURNIN=BURNIN, DECAY_RATE=DECAY_RATE, T_INIT=T_INIT, T_FINAL=T_FINAL,
                     PROB_ADJ=PROB_ADJ, MAKEPLOT=MAKEPLOT, dG_END_LIMIT=dG_END_LIMIT,
                     dG_MID_LIMIT=dG_MID_LIMIT, dG_BAD_LIMIT=dG_BAD_LIMIT)
    CPUs = THREADS if THREADS else 1
    if CPUs>N_RUNS:
        CPUs=N_RUNS
    with Pool(CPUs) as pool:
        # calc number of runs per thread
        chunksize = ceil(N_RUNS / CPUs)
        results = pool.map_async(worker, #function
                                 list(range(1,N_RUNS+1)), #iterable
                                 chunksize=chunksize #max tasks per processor
                                 ).get()
        for row in results:
            loads.append(row)
    
    # export dimer loads per run
    with open(OUTPATH+'_RunSummary.csv', 'w') as file:
        writer = csv.writer(file)
        for row in loads:
            writer.writerow(row)
        
    # remove panel assessment files
    assess_files = glob.glob(OUTPATH+"_Run*_primers_*")
    for f in assess_files:
        os.remove(f)

    # move outputs into separate folders for clarity
    OUTDIR=os.path.dirname(OUTPATH) or "."
    #moveAllFiles(outdir+"/*_SAprimers.csv", os.path.join(newoutdir, "Checkpoint_Primers"))
    #moveAllFiles(outdir+"/*_SAdimers.csv", os.path.join(newoutdir, "Checkpoint_Primers"))
    moveAllFiles(OUTDIR+"/*_primers.csv", os.path.join(OUTDIR, "Final_Primers"))
    moveAllFiles(OUTDIR+"/*_dimers.csv", os.path.join(OUTDIR, "Final_Dimers"))
    moveAllFiles(OUTDIR+"/*_costsTrace.csv", os.path.join(OUTDIR, "Trace_Dimer_Load"))
    moveAllFiles(OUTDIR+"/*_DimerLoad.png", os.path.join(OUTDIR, "Plots_Dimer_Load"))
    moveAllFiles(OUTDIR+"/*log", os.path.join(OUTDIR, "Logfiles"))



def runOpt(run, PRIMER_FA, DIMER_SUMS, DIMER_TABLE, OUTPATH, N_LOCI, 
           KEEPLIST, deltaG, SEED, VERBOSE, SIMPLE, ITERATIONS, CYCLES,
           BURNIN, DECAY_RATE, T_INIT, T_FINAL, PROB_ADJ, MAKEPLOT,#TIMEOUT=5, 
           dG_MID_LIMIT, dG_END_LIMIT, dG_BAD_LIMIT):
    """
    Function for optimizing panels + assessing results for each run
    (to pass to mp.pool)
    """
    # set timer- this helps to timeout runs with infinite loops in the optimization process
    #TIMEOUT = TIMEOUT*60 # convert to seconds
    #start=time.time()
    #while time.time()-start <= TIMEOUT:
    try:
            OUT = OUTPATH+"_Run"+str(run).zfill(2)
            # run optimization
            optimizeMultiplex(PRIMER_FASTA = PRIMER_FA, 
                              DIMER_SUMS=DIMER_SUMS, 
                              DIMER_TABLE=DIMER_TABLE, 
                              OUTPATH = OUT, 
                              N_LOCI=N_LOCI, 
                              KEEPLIST=KEEPLIST, 
                              deltaG=deltaG, 
                              SEED=SEED, 
                              VERBOSE=VERBOSE,
                              SIMPLE=SIMPLE, 
                              ITERATIONS=ITERATIONS, 
                              CYCLES=CYCLES, 
                              BURNIN=BURNIN, 
                              DECAY_RATE=DECAY_RATE, 
                              T_INIT=T_INIT, 
                              T_FINAL=T_FINAL, 
                              PROB_ADJ=PROB_ADJ, 
                              MAKEPLOT=MAKEPLOT,
                              RNG = run*10)
            # run panel assessment
            cost = assessPanel(OUT+"_primers.csv", 
                               ALL_DIMERS_dG=dG_MID_LIMIT,
                               END_DIMERS_dG=dG_END_LIMIT,
                               BAD_DIMERS_dG=dG_BAD_LIMIT)
            cost.insert(0, "Run"+str(run).zfill(2))
            return cost
            #loads.append(cost)
            #run+=1
    except Exception:
            print("AN ERROR OCCURRED IN RUN "+str(run))
            traceback.print_exc()
            print("---")
            return ["Run"+str(run).zfill(2), "ERROR", None, None, None, None, None, None]            #continue
    #else:
    #    raise TimeoutException("in run "+str(run))
    #continue



def moveAllFiles(filegrep, dest):
    #if os.path.exists(dest):
    #    shutil.rmtree(dest)
    if not os.path.exists(dest):
        os.makedirs(dest, exist_ok=True)
    filelist = glob.glob(filegrep)
    for f in filelist:
        try:
            d = os.path.join(dest, os.path.basename(f))
            if os.path.exists(d): os.remove(d)
            shutil.move(f, dest)
        except Exception:
            print("Failed to move "+f)
            pass


class OptimizationWarning(Exception):
    pass
        

class InputError(Exception):
    pass


### Set up timeout exception behavior
class TimeoutException(Exception):   # Custom exception class
    pass

#def timeout_handler(signum, frame):   # Custom signal handler
#    raise TimeoutException


# Change the behavior of SIGALRM
#signal.signal(signal.SIGALRM, timeout_handler)


def parse_args():
    # initialize argparser
    parser = argparse.ArgumentParser()
    # add required arguments
    parser.add_argument("-r", "--runs", type=int, required=True,
                        help="Number of optimization runs")
    parser.add_argument("-f", "--primer_fasta", type=str, required=True,
                        help="Path to FASTA of candidate primers output by primer-design step")
    parser.add_argument("-d", "--dimer_sums", type=str, required=True,
                        help="Path to dimer totals table output by tabulate-dimers step")
    parser.add_argument("-t", "--dimer_table", type=str, required=True,
                        help="Path to pairwise dimers table output by tabulate-dimers step")
    parser.add_argument("-o", "--outpath", type=str, required=True,
                        help="Prefix for output files, including directory")
    parser.add_argument("-n", "--nloci", type=int, required=True,
                        help="Number of primer pairs in optimized panel")
    # add optional arguments
    parser.add_argument("-k", "--keeplist", type=str, default=None,
                        help="Path to FASTA of keeplist primer sequences")
    parser.add_argument("--seed", type=str, default=None,
                        help="Path to previous optimization output to use as starting point")
    parser.add_argument("-s", "--simple", type=int, default=5000,
                        help="Iterations for simple iterative improvement")
    parser.add_argument("-i", "--iter", type=int, default=1000,
                        help="Iterations for each cycle of adaptive simulated annealing (ASA)")
    parser.add_argument("-c", "--cycles", type=int, default=10,
                        help="Number of adaptive simulated annealing cycles")
    parser.add_argument("-b", "--burnin", type=int, default=200,
                        help="Iterations for sampling cost space to set ASA temperature schedule")
    parser.add_argument("--decay-rate", type=float, default=0.95,
                        help="ASA temperature decay rate")
    parser.add_argument("--t-init", type=float, default=None,
                        help="Initial temperature for fixed-schedule simulated annealing")
    parser.add_argument("--t-final", type=float, default=0.01,
                        help="End temperature for fixed-schedule simulated annealing")
    parser.add_argument("--prob-adj", type=float, default=2,
                        help="Adjustment multiplier for ASA acceptance probabilities")
    #parser.add_argument("--timeout", type=float, default=5,
    #                    help="Maximum allowed time (minutes) per swap")
    # add flags
    parser.add_argument("-g", "--deltaG", action="store_true",
                        help="Use deltaG optimization algorithm")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Print out all steps")
    parser.add_argument("-m", "--makeplot", action="store_true",
                        help="Make temperature schedule plots for each run")
    # add sub-module arguments
    parser.add_argument("--threads", type=int, default=None,
                        help="Number of processors for multi-threading")
    parser.add_argument("--dg_end_limit", type=float, default=-4,
                        help="DeltaG threshold (kcal/mol) for 3' end dimers used in panel assessment")
    parser.add_argument("--dg_mid_limit", type=float, default=-8,
                        help="DeltaG threshold (kcal/mol) for non-end dimers used in panel assessment")
    parser.add_argument("--dg_bad_limit", type=float, default=-10,
                        help="DeltaG threshold (kcal/mol) for 'bad' dimers used in panel assessment")
    return parser.parse_args()



def cli():
    # parse command line arguments
    args = parse_args()
    # run main
    main(N_RUNS = args.runs, 
         PRIMER_FA = args.primer_fasta, 
         DIMER_SUMS = args.dimer_sums, 
         DIMER_TABLE = args.dimer_table, 
         OUTPATH = args.outpath, 
         N_LOCI = args.nloci, 
         deltaG = args.deltaG,
         KEEPLIST = args.keeplist, 
         #TIMEOUT = args.timeout, 
         VERBOSE = args.verbose, 
         SEED = args.seed,
         SIMPLE = args.simple, 
         ITERATIONS = args.iter, 
         CYCLES = args.cycles,
         BURNIN = args.burnin, 
         DECAY_RATE = args.decay_rate,
         T_INIT = args.t_init, 
         T_FINAL = args.t_final, 
         PROB_ADJ = args.prob_adj,
         MAKEPLOT = args.makeplot,
         THREADS = args.threads,
         dG_END_LIMIT = args.dg_end_limit,
         dG_MID_LIMIT = args.dg_mid_limit, 
         dG_BAD_LIMIT = args.dg_bad_limit)



if __name__=="__main__":
    cli()
