#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: MULTIPLEX WORMHOLE (Python wrapper)
Purpose: multiplex_wormhole optimizes primer design for multiplex amplicon sequencing 
    by minimizing predicted pairwise dimers. The target audience is for SNP panel
    development, however the process is transferable to any application where multiple
    distinct amplicons are being targeted for multiplex PCR.

Created on Mon Aug 18 22:45:01 2025
@author: maggiehallerud
"""
# load dependencies and modules
import sys
import importlib
import logging
import traceback
from datetime import datetime
import os
import shutil
import pandas as pd
import subprocess
import argparse
import json


# load multiplex wormhole functions
# load within module
try:
    from .batch_primer3_design import main as primer3BatchDesign
    from .tabulate_dimers import main as tabulateDimers
    from .multiple_run_optimization import main as multipleOptimizations
    from .helpers.CSVtoFasta import main as CSVtoFASTA
    from .helpers.logging_setup import setup_logging
    from .helpers._setup_mfeprimer import main as setup_mfeprimer
    from . import plot_ASA_temps as plotASAtemps
    from .optimize_multiplex import LoadPrimers
# load using standalone script
except ImportError:
    sys.path.append(os.path.dirname(__file__))
    from batch_primer3_design import main as primer3BatchDesign
    from tabulate_dimers import main as tabulateDimers
    from multiple_run_optimization import main as multipleOptimizations
    from helpers.CSVtoFasta import main as CSVtoFASTA
    from helpers.logging_setup import setup_logging
    from helpers._setup_mfeprimer import main as setup_mfeprimer
    from optimize_multiplex import LoadPrimers
    plotASAtemps = importlib.import_module("plot_ASA_temps")



def main(TEMPLATES, N_LOCI, OUTDIR, PREFIX=None, KEEPLIST_FA=None, N_RUNS=10, 
         ITERATIONS=1000, CYCLES=10, SIMPLE=5000, deltaG=False, VERBOSE=False, 
         # primer design options                      
         Tm_LIMIT=45, dG_HAIRPINS=-2, dG_END_LIMIT=-4,  dG_MID_LIMIT=-8, 
         ENABLE_BROAD=False, PRIMER3_SETTINGS=None, 
         # MFEprimer options (also uses dG_END_LIMIT & dG_MID_LIMIT from above)
         THREADS=None,
         # optimization / plotting options
         BURNIN=200, DECAY_RATE=0.95, T_INIT=None, T_FINAL=0.01, PROB_ADJ=2,
         # assessment limits
         dG_BAD_LIMIT=-10):
    """
    ----------
    TEMPLATES : CSV containing DNA sequences, IDs, and target sin "startBP,lengthBP" format. [filepath]
    N_LOCI : Target number of amplicons in multiplex PCR. [int]
    OUTDIR : Directory where all outputs will be stored. [Filepath]    
    PREFIX : Prefix to use for optimization output files. [Default: timestamp]
    KEEPLIST_FA : FASTA of primers that MUST be included in final multiplex. [Default: None)
    N_RUNS : Number of optimization runs. [Default: 10]
    ITERATIONS : Iterations per simulated annealing cycle. [Default: 1000]
    CYCLES : Number of simulated annealing cycles to run. [Default: 10]
    SIMPLE : Iterations to run simple iterative improvement optimization. [Default: 5000]
    deltaG : True (deltaG optimization) / False (standard optimization). [Default: False]
    VERBOSE : Print updates as function runs? (Default: False)
    -------
    Returns
    1. Designs & filters primers for each candidate sequences in TEMPLATES
    2. Pairwise dimers predicted between all filtered primer pairs
    3. Optimized multiplex primer set for target N_LOCI
    """    
    # set up MFEprimer
    MFEprimer_PATH = setup_mfeprimer()
    
    # check all inputs
    if not os.path.exists(MFEprimer_PATH):
        raise Exception("MFEprimer_PATH not found! Path provided on lines 37-38: "+MFEprimer_PATH)
    if not os.path.exists(TEMPLATES):
        raise InputError("TEMPLATES file not found!")
    if KEEPLIST_FA is not None:
        if not os.path.exists(KEEPLIST_FA):
            raise InputError("KEEPLIST_FA file not found!")
    # set suffix to current datetime if not given
    if PREFIX is None:
        PREFIX = str(datetime.now()).replace(" ","_").replace(".","_").replace(":","")
    
    # check format of KEEPLIST primerIDs
    try:
        test = LoadPrimers(KEEPLIST_FA)
    except Exception:
        raise InputError("PrimerIDs in the KEEPLIST_FA are not in the proper format! "\
                         "Reformat as <locus>.<#>.<DIR>, e.g., MACA01.0.FWD & MACA01.0.REV")
    
    # initialize logging
    mwlogger = setup_logging("multiplex_wormhole_"+PREFIX+".log", VERBOSE, "mw_main")
    # log start time & inputs
    mwlogger.info("START TIME: %s", datetime.now().strftime('%m/%d/%Y %I:%M:%S %p'))
    mwlogger.info("")
    mwlogger.info("multiplex_wormhole inputs: ")
    mwlogger.info("     TEMPLATES: %s", TEMPLATES)
    mwlogger.info("     N_LOCI: %s", N_LOCI)
    mwlogger.info("     OUTDIR: %s", OUTDIR)
    mwlogger.info("     PREFIX: %s", PREFIX)
    mwlogger.info("     KEEPLIST_FA: %s", KEEPLIST_FA)
    mwlogger.info("     N_RUNS: %s", N_RUNS)
    mwlogger.info("     SIMPLE: %s", SIMPLE)
    mwlogger.info("     ITERATIONS: %s", ITERATIONS)
    mwlogger.info("     deltaG: %s", deltaG)
    mwlogger.info("     VERBOSE: %s", VERBOSE)
    mwlogger.info("     THREADS: %s", THREADS)
    mwlogger.info("     Tm_LIMIT: %s", Tm_LIMIT)
    mwlogger.info("     dG_MID_LIMIT: %s", dG_MID_LIMIT)
    mwlogger.info("     dG_END_LIMIT: %s", dG_END_LIMIT)
    mwlogger.info("     dG_BAD_LIMIT: %s", dG_BAD_LIMIT)
    mwlogger.info("     ENABLE_BROAD: %s", ENABLE_BROAD)
    if PRIMER3_SETTINGS is not None:
        mwlogger.info("     PRIMER3_SETTINGS: (see primerdesign logfile)")
    mwlogger.info("     BURNING: %s", BURNIN)
    mwlogger.info("     DECAY_RATE: %s", DECAY_RATE)
    mwlogger.info("     PROB_ADJ: %s", PROB_ADJ)
    mwlogger.info("     T_INIT: %s", T_INIT)
    mwlogger.info("     T_FINAL: %s", T_FINAL)
    
    
    # set input types
    N_LOCI = int(N_LOCI)
    N_RUNS = int(N_RUNS)
    ITERATIONS = int(ITERATIONS)
    SIMPLE = int(SIMPLE)
    
    ## Step 0: Set up output directory structure & copy inputs to it
    mwlogger.info("-----SETTING UP OUTPUT DIRECTORY STRUCTURE------")
    # set up folder structure
    os.makedirs(OUTDIR, exist_ok=True)
    INPUTDIR = os.path.join(OUTDIR, '0_Inputs')
    os.makedirs(INPUTDIR, exist_ok=True)
    OUTDIR1 = os.path.join(OUTDIR, "1_PrimerDesign")
    os.makedirs(OUTDIR1, exist_ok=True)
    OUTDIR2 = os.path.join(OUTDIR, "2_PredictedDimers")
    os.makedirs(OUTDIR2, exist_ok=True)
    OUTDIR3 = os.path.join(OUTDIR, "3_OptimizedMultiplexes")
    os.makedirs(OUTDIR3, exist_ok=True)
    # copy inputs here
    shutil.copy2(TEMPLATES, INPUTDIR)
    if KEEPLIST_FA is not None:
        shutil.copy2(KEEPLIST_FA, INPUTDIR)
    
    ## Step 1: batch design of primers (with intra-pair hairpin and dimer filtering)
    mwlogger.info("")
    mwlogger.info("-----BATCH DESIGNING PRIMERS------")
    mwlogger.info("Logging to: %s", os.path.join(OUTDIR1, "FilteredPrimers.log"))
    # NOTE: This script includes all of the primer design settings!
    # These can be adjusted directly in the primer3_batch_design.py script
    # or by providing a dictionary to SETTINGS.
    # details on settings: https://htmlpreview.github.io/?https://github.com/primer3-org/primer3/blob/v2.6.1/src/primer3_manual.htm#globalTags 
    primer3BatchDesign(TEMPLATES, 
                       os.path.join(OUTDIR1, "FilteredPrimers"),
                       Tm_LIMIT=Tm_LIMIT, #upper limit for dimer melting temp
                       dG_HAIRPINS=dG_HAIRPINS,  #lower limit for hairpin deltaG
                       dG_END_LIMIT=dG_END_LIMIT, #lower limit for 3' end dimer deltaG
                       dG_MID_LIMIT=dG_MID_LIMIT, #lower limit for deltaG of all other dimers
                       KEEPLIST=KEEPLIST_FA, #keeplist FASTA
                       ENABLE_BROAD=False, #try broader settings if primer design fails?
                       SETTINGS=None) #primer3 settings (in dictionary format)
    # Outputs are found in the 1_PrimerDesign folder. 
    
    
    
    ## Step 2: Predict primer dimers using MFEprimer
    # set input based on whether KEEPLIST is provided or not:
    if KEEPLIST_FA is None:
        INPUT = os.path.join(OUTDIR1, "FilteredPrimers.fa")
    else:
        INPUT = os.path.join(OUTDIR1, "FilteredPrimers_plusKeeplist.fa")
    mwlogger.info("")
    mwlogger.info("-----RUNNING DIMER PREDICTION WITH MFEPRIMER-----")
    mwlogger.info("Note: This is the slowest step! Especially with large candidate sets. Use --threads option for multiprocessing.")
   # NOTE: Originally, primers were checked via the PrimerSuite PrimerDimer function (http://www.primer-dimer.com/)
    # PrimerSuite PrimerDimerReport files can be converted to the necessary table/sum files using scripts/translate_primerSuite_report.R
    # I decided to transition to MFEprimer because primer-dimer.com returned an unreasonable number of dimers
    # set output paths
    ALL_DIMERS=os.path.join(OUTDIR2, 'MFEprimerDimers.txt')
    END_DIMERS=os.path.join(OUTDIR2, 'MFEprimerDimers_ends.txt')
    # MFEprimer parameters:
    # -i = input FASTA of primer sequences 
    # -o = output file
    # -d = maximum deltaG threshold to consider dimers (kcal/mol)
    # -s = minimum score threshold to consider dimers(scores are calculated with +1 for each match and -1 for each mismatch (not including Ns)
    # -m = max allowed mismatches per dimer
    # -p = only output dimers with 3' end bind
    # --diva = concentration of divalent cations (mM)
    # --mono = concentration of monovalent cations (mM)
    # --dntp = concentration of dNTPs (mM)
    # --oligo = concentration of annealing oligos (nM)
    if not os.path.exists(INPUT):
        raise Exception(INPUT+" not found- did primer3BatchDesign step fail?")
    try:
        CPU = THREADS if THREADS else 2
        subprocess.call(MFEprimer_PATH+" dimer -i "+INPUT+" -o "+ALL_DIMERS+" -d "+str(dG_MID_LIMIT)+" -s 3 -m 50 --diva 3.8 --mono 50 --dntp 0.25 --oligo 50 --cpu "+str(CPU),
                        shell=True)
        subprocess.call(MFEprimer_PATH+" dimer -i "+INPUT+" -o "+END_DIMERS+" -d "+str(dG_END_LIMIT)+" -s 3 -m 70 --diva 3.8 --mono 50 --dntp 0.25 --oligo 50 -p --cpu "+str(CPU),
                        shell=True)
    except Exception:
        print("MFEprimer dimer step failed! Error message & traceback:")
        print(traceback.format_exc())
    
    
    
    ## Step 3: Convert MFEprimer dimer report to table formats
    ## NOTE: This is the most computationally intensive step. 
    ## It will run substantially faster if you leave the 4th argument blank 
    ## (which means pairwise interactions between individual primers won't be calculated)
    mwlogger.info("")
    mwlogger.info("-----TABULATING PREDICTED DIMERS-----")
    mwlogger.info("Logging to: %s", os.path.join(OUTDIR2, 'PrimerPairInteractions.log'))
    tabulateDimers(ALL_DIMERS, 
                   END_DIMERS, 
                   os.path.join(OUTDIR2, 'PrimerPairInteractions'), 
                   "False",#os.path.join(OUTDIR3, 'RawPrimerInteractions'))#specify this parameter if you care about per-primer dimers (Rather than just sums per primer pair)
                   deltaG)
    # Outputs are found under 3_PredictedDimers/PrimerPairInteractions*
    if deltaG:
        DIMER_TOTS = os.path.join(OUTDIR2, 'PrimerPairInteractions_mean.csv')
    else:
        DIMER_TOTS = os.path.join(OUTDIR2, 'PrimerPairInteractions_sum.csv')
    
    
    ## Step 4 (Optional): Explore temperature space for simulated annealing
    ## There are two ways to run this script: one calculates temperatures and dimer loads based on the problem at hand, 
    ## the other uses pre-specified temperatures and dimer loads.
    ## I recommend first running using files from the problem, then using the values observed in the outputs to explore 
    ## parameters around the defaults.
    mwlogger.info("")
    mwlogger.info("-----PLOTTING ADAPTIVE SIMULATED ANNEALING TEMPERATURE SCHEDULE------")
    mwlogger.info("Logging to: %s", os.path.join(OUTDIR3, 'ASAplots.log'))
    plotASAtemps.main(OUTPATH=os.path.join(OUTDIR3, 'ASAplots'),
                      PRIMER_FASTA=os.path.join(OUTDIR1, 'FilteredPrimers.fa'), 
                      DIMER_SUMS=DIMER_TOTS,
                      DIMER_TABLE=os.path.join(OUTDIR2, 'PrimerPairInteractions_wide.csv'), 
                      N_LOCI=N_LOCI, #number of target loci in panel
                      KEEPLIST=KEEPLIST_FA, 
                      SEED=None, #this would be an output from optimizeMultiplex
                      BURNIN=BURNIN,#number iterations with dimer loads used to sample cost space
                      DECAY_RATE=DECAY_RATE, #temp decay rate
                      PROB_ADJ=PROB_ADJ,#decay rate of acceptance probability
                      deltaG=deltaG,
                      T_INIT=T_INIT,
                      T_FINAL=T_FINAL)
    # Alternative implementation where dimers and T_INIT/T_FINAL are specified:
    # plotSAtemps(OUTPATH=os.path.join(OUTDIR4, "TestingASAparams"),
    #             # dimer counts to plot and calculate temps from (if not set)
    #             MIN_DIMER=1,
    #             MAX_DIMER=5, #update this value based on the max observed in the default plot!
    #             # parameter determining temperature decay in negative exponential
    #             DECAY_RATE=0.95, #default=0.95
    #             # initial temperature to start from - higher=more risk accepted
    #             T_INIT=2, #default=2
    #             # final temperature to stop at - 1=no risk, higher=more risk accepted
    #             T_FINAL=0.01, 
    #             #proportion of max dimer load considered when setting temperature schedule
    #             #closer to 0 = accepts fewer errors
    #             # adjustment for dimer acceptance probabilities- 1=no adjustment, higher values=lower dimer acceptance
    #             PROB_ADJ=1)#default=2
    
    
    ## Step 5: Design a set of multiplex primers by minimizing predicted dimer formation
    # N_LOCI here is the number of loci you want in the final panel (including keeplist loci) 
    # N_RUNS is number of runs of optimization process (1 for simple problems, increase for complex problems)
    mwlogger.info("")
    mwlogger.info("-----STARTING OPTIMIZATION FOR MULTIPLEX PRIMER SET-----")
    mwlogger.info("Logged to %s",  os.path.join(OUTDIR3, PREFIX)+".log")
    multipleOptimizations(N_RUNS = N_RUNS, 
                          PRIMER_FA = os.path.join(OUTDIR1, "FilteredPrimers.fa"),
                          DIMER_SUMS = DIMER_TOTS,
                          DIMER_TABLE = os.path.join(OUTDIR2, 'PrimerPairInteractions_wide.csv'), 
                          OUTPATH = os.path.join(OUTDIR3, PREFIX),
                          N_LOCI = N_LOCI, 
                          deltaG = deltaG, #True: deltaG optimization, False=standard optimization
                          KEEPLIST = KEEPLIST_FA, 
                          TIMEOUT = 360,#time allowed per run- runs 30 minutes will break
                          VERBOSE=VERBOSE,#set to true to print dimers at each change
                          SIMPLE=SIMPLE, # iterations for simple iterative improvement optimization (default=5000)
                          ITERATIONS=ITERATIONS, # iterations per simulated annealing optimization cycle (default=1000) 
                          CYCLES=CYCLES, #simulated annealing cycles to run (default=10)
                          BURNIN=BURNIN, # iterations for sampling dimer cost space to adaptively set SA temps (default=100)
                          DECAY_RATE=DECAY_RATE, # temperature decay parameter for SA temps (default=0.95)
                              # closer to 1 - least conservative, explores more cost space at higher risk
                              # closer to 0 - most conservative, explores less cost space at lower risk
                              # recommendations: 0.90-0.98, higher with fewer iterations
                          T_INIT=T_INIT, # starting temp for fixed SA schedule (default=None, i.e., adaptively set based on costs observed in BURNIN)
                          T_FINAL=T_FINAL, # ending temp for fixed SA schedule (default=0.1)
                              # temperatures=0 is equivalent to simple iterative improvement, while 
                              # higher temperatures explore more of the cost space at higher risk of accepting dimers
                              # recommended initial fixed schedule is T_INIT~2 and T_FINAL=0.1
                          PROB_ADJ=PROB_ADJ,# adjusts dimer acceptance probabilities (default=2)
                              # increase if too many dimers are being accepted during simulated annealing, 
                              # decrease if local optima are not being overcome
                          SEED=None,#primer set from previous optimization run to start with, in CSV format
                          # inputs to mw.assessPanel
                          dG_END_LIMIT=dG_END_LIMIT,
                          dG_MID_LIMIT=dG_MID_LIMIT,
                          dG_BAD_LIMIT=dG_BAD_LIMIT)
    #OUTPUT: MAF30_150loci_RunSummary.csv
    
    
    ## STEP 7: Convert selected primer set to FASTA format for additional screening
    runs = pd.read_csv(os.path.join(OUTDIR3, PREFIX+"_RunSummary.csv"))
    if len(runs)>1:
        mwlogger.info("")
        mwlogger.info("-----CONVERTING BEST PRIMER SETS INTO FASTAs FOR ADDITIONAL SCREENING-----")
        run = runs['Run']
        run = [str(run[x]).zfill(2) for x in range(len(run))]
        dimers = runs['TotalDimers']
        mwlogger.info("The BEST multiplex had %s total predicted dimers.", str(min(dimers)))
        for i in range(len(runs)):
            mwlogger.info("Run %s had %s dimers.", str(run[i]), str(dimers[i]))
            if dimers[i]==min(dimers):
                mwlogger.info(".....Converting to FASTA for additional screening")
                CSVtoFASTA(IN_CSV = os.path.join(OUTDIR3,"Final_Primers", PREFIX+"_"+run[i]+"_primers.csv"), 
                           OUT_FA = os.path.join(OUTDIR3,PREFIX+run[i]+"_primers.fasta"))
    else: 
        mwlogger.warning("RunSummary not found- optimization step may have failed.")
    
    # finish logging
    # end logging
    mwlogger.info("END TIME: %s", datetime.now().strftime('%m/%d/%Y %I:%M:%S %p'))
    print("")
    print("LOGGED TO: multiplex_wormhole_"+PREFIX+".log")
    logging.shutdown()
    


class InputError(Exception):
    pass



def parse_args():
    # initialize argparser
    parser = argparse.ArgumentParser()
    # add required arguments
    parser.add_argument("-t", "--templates", type=str, required=True, 
                        help="Filepath to CSV containing SEQUENCE_ID, SEQUENCE_TEMPLATE, and SEQUENCE_TARGET fields")
    parser.add_argument("-n", "--nloci", type=int, required=True,
                        help="Number of primer pairs in optimized multiplex")
    parser.add_argument("-o", "--outdir", type=str, required=True,
                        help="Name of output directory to create & store outputs")
    # add optional arguments
    parser.add_argument("-p", "--prefix", type=str, default=None,
                        help="Prefix used to name output files")
    parser.add_argument("-k", "--keeplist", type=str, default=None,
                        help="Filepath to FASTA containing keeplist primer sequences")
    parser.add_argument("-r", "--runs", type=int, default=10,
                        help="Number of optimization runs")
    parser.add_argument("-s", "--simple", type=int, default=5000,
                        help="Number of iterations in simple iterative improvement")
    parser.add_argument("-i", "--iter", type=int, default=1000,
                        help="Number of iterations in each adaptive simulated annealing (ASA) cycle")
    parser.add_argument("-c", "--cycles", type=int, default=10,
                        help="Number of adaptive simulated annealing cycles")
    parser.add_argument("-d", "--deltaG", action="store_true",
                        help="Use deltaG optimization algorithm")#type=str, default=False)
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose logging")#type=str, default=False)
    # add arguments from sub-modules
    parser.add_argument("--tm-limit", type=float, default=45,
                        help="Min melting temperature of dimers (Celsius)")
    parser.add_argument("--dg-hairpins", type=float, default=-2,
                        help="DeltaG threshold for primer hairpins (kcal/mol)")
    parser.add_argument("--dg-end-limit", type=float, default=-4,
                        help="DeltaG threshold for 3' end dimers (kcal/mol)")
    parser.add_argument("--dg-mid-limit", type=float, default=-8,
                        help="DeltaG threshold for non-end dimers (kcal/mol)")
    parser.add_argument("--dg-bad-limit", type=float, default=-10,
                        help="DeltaG threshold for counting 'bad' dimers (kcal/mol)")
    parser.add_argument("--enable-broad", action="store_true",
                        help="Enable broader settings during primer design")
    parser.add_argument("--primer3-settings", type=json.loads, default=None,
                        help="Dictionary of primer3 settings for primer design")
    parser.add_argument("--threads", type=int, default=None,
                        help="# Processors for multi-threading")
    parser.add_argument("--burnin", type=int, default=200,
                        help="Sampling iterations for setting ASA temperatures")
    parser.add_argument("--decay-rate", type=float, default=0.95,
                        help="ASA temperature decay rate")
    parser.add_argument("--t-init", type=float, default=None,
                        help="Initial temperature for fixed-schedule simulated annealing")
    parser.add_argument("--t-final", type=float, default=0.01,
                        help="Final temperature for fixed-schedule simulated annealing")
    parser.add_argument("--prob-adj", type=float, default=2,
                        help="Adjustment for acceptance probabilities in ASA algorithm")
    return parser.parse_args()



def cli():
    # parse command-line arguments
    args = parse_args()
    # run multiplex wormhole module
    main(TEMPLATES = args.templates,
         N_LOCI = args.nloci,
         OUTDIR = args.outdir,
         PREFIX = args.prefix,
         KEEPLIST_FA = args.keeplist,
         N_RUNS = args.runs,
         ITERATIONS = args.iter,
         CYCLES = args.cycles,
         SIMPLE = args.simple,
         deltaG = args.deltaG,
         VERBOSE = args.verbose,
         Tm_LIMIT=args.tm_limit,
         dG_HAIRPINS=args.dg_hairpins,
         dG_END_LIMIT=args.dg_end_limit,
         dG_MID_LIMIT=args.dg_mid_limit,
         dG_BAD_LIMIT=args.dg_bad_limit,
         ENABLE_BROAD=args.enable_broad,
         PRIMER3_SETTINGS=args.primer3_settings,
         THREADS=args.threads,
         BURNIN=args.burnin, 
         DECAY_RATE=args.decay_rate, 
         T_INIT=args.t_init, 
         T_FINAL=args.t_final, 
         PROB_ADJ=args.prob_adj)



if __name__ == "__main__":
    cli()
