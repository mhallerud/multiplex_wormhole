#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: MULTIPLEX WORMHOLE (Python wrapper)
Created on Mon Aug 18 22:45:01 2025
@author: maggiehallerud

Purpose: multiplex_wormhole optimizes primer design for multiplex amplicon sequencing 
    by minimizing predicted pairwise dimers. The target audience is for SNP panel
    development, however the process is transferable to any application where multiple
    distinct amplicons are being targeted for multiplex PCR.
    
IMPORTANT NOTE: Before running this script, make sure that you have installed the 
dependencies primer3 (available at https://github.com/primer3-org/primer3/releases)
and MFEprimer (available at https://www.mfeprimer.com/mfeprimer-3.1/#2-command-line-version)

Input preparation:
    TEMPLATES : A CSV containing loci being targeted for multiplex amplicon sequencing
        The TEMPLATES csv should be in openprimer format:
        - 3 columns: locus ID, DNA sequence, and target position (start bp, length following primer3 format)
        - each target locus should have a separate entry in the CSV
        - locus IDs MUST be unique for multiplex_wormhole to work correctly!
        The script create_in_templates.R will create the template CSV by taking an input VCF file 
        containing target SNPs and a FASTA file containing the sequences associated with these SNPs 
        (with FASTA seq IDs matching the CHROM field in the VCF)
        
    KEEPLIST_FA : A FASTA formatted file containing adapter-ligated primer sequences from a current multiplex
        assay that you are looking to add to. 
        Forward primers must be designated with ".FWD" or ".FW" as a suffix and reverse primers with ".REV"
    
    Note that locus IDs in the TEMPLATES file and primer IDs in the KEEPLIST_FA must be unique
    for multiplex_wormhole to function properly. IDs may not contain periods.
"""
# load dependencies and modules
import sys
import importlib
import glob
import os
import shutil
import pandas as pd
import subprocess
import datetime
import argparse

# load multiplex wormhole functions
sys.path.append(os.path.dirname(__file__))
from batch_primer3_design import main as primer3BatchDesign
from tabulate_dimers import main as tabulateDimers
from multiple_run_optimization import main as multipleOptimizations
from CSVtoFasta import main as CSVtoFASTA
plotASAtemps = importlib.import_module("plot_ASA_temps")

## FINE PATH TO BINARY DEPENDENCIES
## NO SPACES ALLOWED IN PATHS- OTHERWISE CALLING FUNCTIONS WILL BREAK!
MFEprimer_PATH = glob.glob(os.path.dirname(__file__)+"/*mfeprimer*")[0]



def main(TEMPLATES, N_LOCI, OUTDIR, PREFIX=None, KEEPLIST_FA=None, N_RUNS=10, \
         ITERATIONS=10000, SIMPLE=5000, deltaG=False, VERBOSE=False):
    """
    NOTE: DEPENDENCY PATHS MUST BE SET IN THIS SCRIPT FOR IT TO RUN
    Parameters
    ----------
    TEMPLATES : CSV filepath
        Contains ID, DNA sequence, and target in "startBP,lengthBP" format
    N_LOCI : integer
        Target number of amplicons in multiplex PCR
    OUTDIR : Filepath
        Directory where all outputs will be stored
    PREFIX : String
        Prefix to use for optimization output files
    KEEPLIST_FA : FASTA filepath, optional (default: None)
        Primers that MUST be included in final multiplex set.
    N_RUNS : integer (default: 10)
        Number of optimization runs
    ITERATIONS : integer (default: 5000)
        Iterations to run simulated annealing optimization
    SIMPLE : integer (default: 2000)
        Iterations to run simple iterative improvement optimization
    deltaG : True (deltaG optimization) / False (standard optimization)
    VERBOSE : print updates as function runs? (Default: False)
    -------
    Returns
    1. Designs primers for each candidate sequences in TEMPLATES
    2. Filtered primer sets in CSV format
    3. Pairwise dimers predicted between all filtered primer pairs
    4. Optimized multiplex primer set for target N_LOCI
    5. Predicted dimer loads for primers included in optimized multiplex
    6. Plots of simulated annealing temperature schedule and trace of dimer load
    """
    
    # check for MFEprimer path
    if not os.path.exists(MFEprimer_PATH):
        raise Exception("MFEprimer_PATH not found! Path provided on lines 37-38: "+MFEprimer_PATH)
    
    # set input types
    N_LOCI = int(N_LOCI)
    N_RUNS = int(N_RUNS)
    ITERATIONS = int(ITERATIONS)
    SIMPLE = int(SIMPLE)
    
    ## Step 0: Set up output directory structure & copy inputs to it
    print("-----SETTING UP OUTPUT DIRECTORY STRUCTURE------")
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
    
    # set suffix to current datetime if not given
    if PREFIX is None:
        PREFIX = str(datetime.datetime.now()).replace(" ","_").replace(".","_")
        
    ## Step 1: batch design of primers (with intra-pair hairpin and dimer filtering)
    print("")
    print("-----BATCH DESIGNING PRIMERS------")
    # NOTE: This script includes all of the primer design settings!
    # These can be adjusted directly in the primer3_batch_design.py script
    # or by providing a dictionary to SETTINGS.
    # details on settings: https://htmlpreview.github.io/?https://github.com/primer3-org/primer3/blob/v2.6.1/src/primer3_manual.htm#globalTags 
    primer3BatchDesign(TEMPLATES, 
                       os.path.join(OUTDIR1, "FilteredPrimers"),
                       Tm_LIMIT=45, #upper limit for dimer melting temp
                       dG_HAIRPINS=-2,  #lower limit for hairpin deltaG
                       dG_END_LIMIT=-4, #lower limit for 3' end dimer deltaG
                       dG_MID_LIMIT=-8, #lower limit for deltaG of all other dimers
                       KEEPLIST=None, #keeplist FASTA
                       ENABLE_BROAD=False, #try broader settings if primer design fails?
                       SETTINGS=None) #primer3 settings (in dictionary format)
    # Outputs are found in the 1_PrimerDesign folder. 
    
    
    
    ## Step 2: Predict primer dimers using MFEprimer
    # set input based on whether KEEPLIST is provided or not:
    if KEEPLIST_FA is None:
        INPUT = os.path.join(OUTDIR1, "FilteredPrimers.fa")
    else:
        INPUT = os.path.join(OUTDIR1, "FilteredPrimers_plusKeeplist.fa")
    print("")
    print("-----RUNNING DIMER PREDICTION WITH MFEPRIMER-----")
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
    subprocess.call(MFEprimer_PATH+" dimer -i "+INPUT+" -o "+ALL_DIMERS+" -d -8 -s 3 -m 50 --diva 3.8 --mono 50 --dntp 0.25 --oligo 50",
                     shell=True)
    subprocess.call(MFEprimer_PATH+" dimer -i "+INPUT+" -o "+END_DIMERS+" -d -3 -s 3 -m 70 --diva 3.8 --mono 50 --dntp 0.25 --oligo 50 -p",
                    shell=True) 
    
    
    
    ## Step 3: Convert MFEprimer dimer report to table formats
    ## NOTE: This is the most computationally intensive step. 
    ## It will run substantially faster if you leave the 4th argument blank 
    ## (which means pairwise interactions between individual primers won't be calculated)
    print("")
    print("-----TABULATING PREDICTED DIMERS-----")
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
    print("")
    print("-----PLOTTING ADAPTIVE SIMULATED ANNEALING TEMPERATURE SCHEDULE------")
    plotASAtemps.main(OUTPATH=os.path.join(OUTDIR3, 'ASAplots'),
                      PRIMER_FASTA=os.path.join(OUTDIR1, 'FilteredPrimers.fa'), 
                      DIMER_SUMS=DIMER_TOTS,
                      DIMER_TABLE=os.path.join(OUTDIR2, 'PrimerPairInteractions_wide.csv'), 
                      N_LOCI=N_LOCI, #number of target loci in panel
                      KEEPLIST=KEEPLIST_FA, 
                      SEED=None, #this would be an output from optimizeMultiplex
                      BURNIN=100,#number iterations with dimer loads used to sample cost space
                      DECAY_RATE=0.95, #temp decay rate
                      DIMER_ADJ=0.1, #adjustment of maximum dimer costs
                      PROB_ADJ=2,#decay rate of acceptance probability
                      deltaG=deltaG)
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
    #             T_FINAL=0, #default=0.1
    #             #proportion of max dimer load considered when setting temperature schedule
    #             #closer to 0 = accepts fewer errors
    #             DIMER_ADJ=0.1, #default=1
    #             # adjustment for dimer acceptance probabilities- 1=no adjustment, higher values=lower dimer acceptance
    #             PROB_ADJ=1)#default=2
    
    
    ## Step 5: Design a set of multiplex primers by minimizing predicted dimer formation
    # N_LOCI here is the number of loci you want in the final panel (including keeplist loci) 
    # N_RUNS is number of runs of optimization process (1 for simple problems, increase for complex problems)
    print("")
    print("-----STARTING OPTIMIZATION FOR MULTIPLEX PRIMER SET-----")
    multipleOptimizations(N_RUNS = N_RUNS, 
                          PRIMER_FA = os.path.join(OUTDIR1, "FilteredPrimers.fa"),
                          DIMER_SUMS = os.path.join(OUTDIR2, 'PrimerPairInteractions_binary_sum.csv'), 
                          DIMER_TABLE = os.path.join(OUTDIR2, 'PrimerPairInteractions_binary_wide.csv'), 
                          OUTPATH = os.path.join(OUTDIR3, PREFIX),
                          N_LOCI = N_LOCI, 
                          deltaG = deltaG, #True: deltaG optimization, False=standard optimization
                          KEEPLIST = KEEPLIST_FA, 
                          TIMEOUT = 360,#time allowed per run- runs 30 minutes will break
                          VERBOSE=VERBOSE,#set to true to print dimers at each change
                          SIMPLE=SIMPLE, # iterations for simple iterative improvement optimization (default=5000)
                          ITERATIONS=ITERATIONS, # iterations for simulated annealing optimization (default=10000) 
                          BURNIN=100, # iterations for sampling dimer cost space to adaptively set SA temps (default=100)
                          DECAY_RATE=0.95, # temperature decay parameter for SA temps (default=0.95)
                              # closer to 1 - least conservative, explores more cost space at higher risk
                              # closer to 0 - most conservative, explores less cost space at lower risk
                              # recommendations: 0.90-0.98, higher with fewer iterations
                          T_INIT=None, # starting temp for fixed SA schedule (default=None, i.e., adaptively set based on costs observed in BURNIN)
                          T_FINAL=0.1, # ending temp for fixed SA schedule (default=0.1)
                              # temperatures=0 is equivalent to simple iterative improvement, while 
                              # higher temperatures explore more of the cost space at higher risk of accepting dimers
                              # recommended initial fixed schedule is T_INIT~2 and T_FINAL=0.1
                          PARTITIONS=1000, # number of times to change temperature schedule (default=1000)
                          DIMER_ADJ=0.1, # proportion of max observed dimer load to consider when setting SA temps (default=0.1)
                              # values closer to 1 will create a temp schedule with higher risk / more cost exploration
                              # values closer to 0 will create a temp schedule with lower risk / less cost exploration
                          PROB_ADJ=2,# adjusts dimer acceptance probabilities (default=2)
                              # increase if too many dimers are being accepted during simulated annealing, 
                              # decrease if local optima are not being overcome
                          SEED=None)#primer set from previous optimization run to start with, in CSV format
    #OUTPUT: MAF30_150loci_RunSummary.csv
    
    
    ## STEP 7: Convert selected primer set to FASTA format for additional screening
    print("")
    print("-----CONVERTING BEST PRIMER SETS INTO FASTAs FOR ADDITIONAL SCREENING-----")
    runs = pd.read_csv(os.path.join(OUTDIR3, PREFIX+"_RunSummary.csv"))
    run = runs['Run']
    run = [str(run[x]).zfill(2) for x in range(len(run))]
    dimers = runs['TotalDimers']
    print("The BEST multiplex had "+str(min(dimers))+ " total predicted dimers.")
    for i in range(len(runs)):
        print("Run "+str(run[i])+" had "+str(dimers[i])+" dimers.")
        if dimers[i]==min(dimers):
            print(".....Converting to FASTA for additional screening")
            CSVtoFASTA(IN_CSV = os.path.join(OUTDIR3,"Final_Primers", PREFIX+"_Run"+str(run[i])+"_primers.csv"), 
                       OUT_FA = os.path.join(OUTDIR3,"Run"+str(run[i])+"_"+PREFIX+"_primers.fasta"))



def usage():
    #print("\nWelcome to multiplex wormhole! Here's how to use the function:\n")
    print("Usage: '+sys.argv[0]+' -t <templatesCSV> -n <#loci> -o <outdir> [-p <prefix> -k <keeplistFA> "+\
          "-r <nruns> -i <iterations> -s <simple> -d <deltaG> -v <verbose> -h <help>]")
    print("\nAdditional documentation can be found at https://github.com/mhallerud/multiplex_wormhole")



def parse_args():
    # initialize argparser
    parser = argparse.ArgumentParser()
    # add required arguments
    parser.add_argument("-t", "--templates", type=str, required=True)
    parser.add_argument("-n", "--nloci", type=int, required=True)
    parser.add_argument("-o", "--outdir", type=str, required=True)
    # add optional arguments
    parser.add_argument("-p", "--prefix", type=str, default="None")
    parser.add_argument("-k", "--keeplist", type=str, default=None)
    parser.add_argument("-r", "--runs", type=int, default=10)
    parser.add_argument("-i", "--iter", type=int, default=10000)
    parser.add_argument("-s", "--simple", type=int, default=5000)
    parser.add_argument("-d", "--deltaG", action="store_true")#type=str, default=False)
    parser.add_argument("-v", "--verbose", action="store_true")#type=str, default=False)

    return parser.parse_args()



if __name__ == "__main__":
    # parse command-line arguments
    args = parse_args()
    # print to standard out
    print("TEMPLATES: "+args.templates)
    print("NLOCI: "+str(args.nloci))
    print("OUTDIR:"+args.outdir)
    print("PREFIX: "+str(args.prefix))
    if args.keeplist is not None:
        print("KEEPLIST: "+args.keeplist)
    else:
        print("KEEPLIST: None")
    print("RUNS: "+str(args.runs))
    print("ITERATIONS: "+str(args.iter))
    print("SIMPLE: "+str(args.simple))
    print("VERBOSE: "+str(args.verbose))
    print("DeltaG: "+str(args.deltaG))
    # run multiplex wormhole module
    main(TEMPLATES = args.templates,
         N_LOCI = args.nloci,
         OUTDIR = args.outdir,
         PREFIX = args.prefix,
         KEEPLIST_FA = args.keeplist,
         N_RUNS = args.runs,
         ITERATIONS = args.iter,
         SIMPLE = args.simple,
         deltaG = args.deltaG,
         VERBOSE = args.verbose)
