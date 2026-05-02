#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: MULTIPLEX PRIMER DESIGN & ASSESSMENT with multiplex_wormhole
Full Documentation: https://mhallerud.github.io/multiplex_wormhole/
Purpose: multiplex_wormhole optimizes primer design for multiplex amplicon sequencing 
    by minimizing predicted pairwise dimers, and assesses existing multiplexes. 
    The target audience is for SNP panel development, however the process is transferable 
    to other targeted sequencing applications.

Dependencies: primer3-py (developed w/ v2.0.0)
              pandas (developed w/ v1.4.4)
              matplotlib (developed w/ v3.5.2)
              MFEprimer v3.2.7
              
Created on Tue Dec 19 19:50:45 2023
@author: maggiehallerud
"""

#-----LOAD MULTIPLEX WORMHOLE FUNCTIONS & DEPENDENCIES------#
import os
import shutil
import subprocess
import multiplex_wormhole as mw

## Alternative: load from locally cloned repository
#MULTIPLEX_WORMHOLE = "/Users/maggiehallerud/Desktop/multiplex_wormhole/" #no spaces!
#import sys
#sys.path.append(MULTIPLEX_WORMHOLE+"/src/multiplex_wormhole")

# load multiplex wormhole functions from cloned repo
#import importlib
#from batch_primer3_design import main as primer3BatchDesign
#from tabulate_dimers import main as tabulateDimers
#from optimize_multiplex import main as optimizeMultiplex
#from multiple_run_optimization import main as multipleOptimizations
#from helpers.CSVtoFasta import main as CSVtoFASTA
#from panel_assessment import main as assessPanel
#plotASAtemps = importlib.import_module("plot_ASA_temps")



#---------------------PANEL ASSESSMENT WORKFLOW ------------------------#
# INPUT: FASTA or CSV (PrimerID, Sequence) file with primers, PrimerIDs following rules
# described below (e.g., "MACA1.FWD" & "MACA1.REV")
mw.assessPanel("Primers.fasta",
               ALL_DIMERS_dG=-8, #deltaG threshold used to predict non-end dimers
               END_DIMERS_dG=-4, #deltaG threshold used to predict 3' end dimers
               BAD_DIMERS_dG=-10) #deltaG threshold for counting "bad" dimers
# OUTPUTS: # primer pairs, total # pairwise dimers, # primer pairs forming dimers
# also dimer files from MFEprimer and dimer tables from tabulateDimers

    
    
##----------PRIMARY WORKFLOW: OPTIMIZED PANEL DESIGN FOR MULTIPLEX PCR-------##
##--------This will run all of the sub-modules under default settings--------##
"""
INPUT PREPARATION:
1) TEMPLATES : A CSV containing information on the DNA template targeted for sequencing.
The following format is required!
    - 3 fields with the names: 
        SEQUENCE_ID : contains names of targets, without punctutation
        SEQUENCE_TEMPLATE : contains DNA template sequence, including target and flanking regions.
        SEQUENCE_TARGET : target region within template, in <start_bp,length> format. 
            for example, a SNP at bp 45 would be 45,1
            a microhaplotype with SNPs at bp 45 & 80 would be 45,35
    - primer pairs will be developed and optimized for each row, so names and sequences must be unique

A helper script (create_in_templates.R) is provided to help convert a VCF + FASTA inputs
into this CSV format. 
            
2) KEEPLIST_FA : A FASTA formatted file containing adapter-ligated primer sequences from a current multiplex
    assay that you are looking to add to. 
    Forward primers must be designated with ".FWD" or ".FW" as a suffix and reverse primers with ".REV"
"""        
## SET MINIMUM INPUTS:
TEMPLATES = "/multiplex_wormhole/examples/Input_Templates.csv"#DNA templates & targets
N_LOCI = 50 # target panel size (# sequences amplified)
OUTDIR = "MW_TEST" # output name
KEEPLIST_FA = None # FASTA containing previously designed primer set
DELTAG = False #set to True if you want to use deltaG optimization algorithm


mw.multiplexWormhole(TEMPLATES = TEMPLATES,
                     N_LOCI = N_LOCI, 
                     OUTDIR = OUTDIR, 
                     PREFIX= "MW_test", 
                     KEEPLIST_FA= KEEPLIST_FA, 
                     N_RUNS=10, # number of full optimization runs
                     ITERATIONS=1000, # iterations per simulated annealing cycle
                     CYCLES=10, # number of simulated annealing cycles
                     SIMPLE=5000, # iterations for simple iterative improvement
                     deltaG=False, # minimize dimer count (False) or mean deltaG of dimers (True)
                     VERBOSE=False)



##-----------STEP-BY-STEP: OPTIMIZED PANEL DESIGN FOR MULTIPLEX PCR----------##
##-----------This allows ultimate flexibility for exploration----------------##
#### STEP 0: SETUP OUTDIR STRUCTURE
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

    
#### STEP 1: BATCH DESIGN PRIMERS
# NOTE: This script includes all of the primer design settings!
# Setting defaults found here: https://mhallerud.github.io/multiplex_wormhole/primer-design
# full settings are printed to the 1_PrimerDesign/{OUTPATH}.log file after running this function
mw.primer3BatchDesign(TEMPLATES, 
                      OUTPATH = os.path.join(OUTDIR1, "FilteredPrimers"), 
                      Tm_LIMIT = 45, #upper threshold for dimer melting temperatures
                      dG_HAIRPINS = -2, #lower threshold for hairpin structure deltaG
                      dG_END_LIMIT = -4, #threshold for 3' end dimer deltaG
                      dG_MID_LIMIT = -8, #threshold for non-end dimer deltaG
                      KEEPLIST = KEEPLIST_FA, 
                      ENABLE_BROAD = False, #use broader settings if no primers for template?
                      SETTINGS = None) #primer3 settings in dictionary {} format, overrides defaults 
## Outputs will be saved to the {OUTDIR}/1_PrimerDesign folder.


#### STEP 2: PREDICT DIMERS USING MFEprimer
# Set input based on whether keeplist is provided or not
if KEEPLIST_FA is None:
    INPUT = os.path.join(OUTDIR1, "FilteredPrimers.fa")
else:
    INPUT = os.path.join(OUTDIR1, "FilteredPrimers_plusKeeplist.fa")

# set output paths
ALL_DIMERS=os.path.join(OUTDIR2, 'MFEprimerDimers.txt')
END_DIMERS=os.path.join(OUTDIR2, 'MFEprimerDimers_ends.txt')

# check mfeprimer path
if not os.path.exists(mw.MFEprimer_PATH):
    mw.MFEprimer_PATH = mw.setup_mfeprimer()
# if this fails, you will need to manually download MFEprimer & set the path
# instructions here: https://mhallerud.github.io/multiplex_wormhole/#installation

# run MFEprimer dimer function for dimer prediction
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
subprocess.call(mw.MFEprimer_PATH+" dimer -i "+INPUT+" -o "+ALL_DIMERS+" -d -8 -s 3 -m 50 --diva 3.8 "+
          "--mono 50 --dntp 0.25 --oligo 50", shell=True)
subprocess.call(mw.MFEprimer_PATH+" dimer -i "+INPUT+" -o "+END_DIMERS+" -d -4 -s 3 -m 70 --diva 3.8 "+
          "--mono 50 --dntp 0.25 --oligo 50 -p", shell=True)
## Outputs will be saved to the {OUTDIR}/2_PredictedDimers folder.
## NOTE: This step will be slow if there are a lot of candidate primer pairs!!


#### STEP 3: CONVERT MFEprimer TEXT OUTPUT TO TABULAR FORMATS
# If deltaG=False (default), outputs will be the count of dimers per pairwise interaction,
# and summarized per primer pair as the sum across all other pairs.
# If deltaG=True, outputs will be the minimum (i.e. worst) deltaG of dimers in each pairwise interaction,
# and summarized per primer pair as the mean across all pairwise interactions.
mw.tabulateDimers(ALL_DIMERS, 
                  END_DIMERS, 
                  OUTPATH = os.path.join(OUTDIR2, 'PrimerPairInteractions'), 
                  OUTPRIMERPATH = "False",#specify a filepath here if you want to output individual primer interactions
                  deltaG = DELTAG)
## Outputs are saved to 2_PredictedDimers/PrimerPairInteractions*
# names for next step will depend on whether you're using the deltaG algorithm or not...
if DELTAG:
    DIMER_TOTS = os.path.join(OUTDIR2, 'PrimerPairInteractions_mean.csv')
else:
    DIMER_TOTS = os.path.join(OUTDIR2, 'PrimerPairInteractions_sum.csv')


#### STEP 4A (OPTIONAL): EXPLORE SIMULATED ANNEALING PARAMETER SETTINGS
## There are two ways to run this script: 
## 1. Calculate temperatures and dimer loads based on the problem at hand.
## 2. Use pre-specified temperatures and dimer loads.
## I recommend trying (1) first, then using this to inform (2) if you want to 
## explore alternatives beyond the defaults.
mw.plotASAtemps(OUTPATH=os.path.join(OUTDIR3, 'TestingASAparams_defaults'),
                PRIMER_FASTA=os.path.join(OUTDIR1, 'FilteredPrimers.fa'), 
                DIMER_SUMS=DIMER_TOTS,
                DIMER_TABLE=os.path.join(OUTDIR2, 'PrimerPairInteractions_wide.csv'), 
                N_LOCI=N_LOCI, #number of target loci in panel
                KEEPLIST=KEEPLIST_FA, 
                SEED=None, #this would be an output from optimizeMultiplex
                BURNIN=100,#number iterations with dimer loads used to sample cost space
                deltaG=DELTAG)
# decay rate closer to 1: 
mw.plotASAtemps(OUTPATH=os.path.join(OUTDIR3, 'TestingASAparams_decayRate98'),
                # dimer counts to plot and calculate temps from (if not set)
                MIN_DIMER=1,
                MAX_DIMER=5, #update based on the max observed in the default plot!
                DECAY_RATE=0.98, # base for geometric functino of temperature decay
                T_INIT=3, #initial A.S.A. temperature: higher=more risk accepted
                T_FINAL=0.01, #final A.S.A. temp: 0=no risk, higher=accepting risk
                #proportion of max dimer load considered when setting temperature schedule
                #1=no adjustment, closer to 0 = less mistake-tolerant
                DIMER_ADJ=0.1,
                PROB_ADJ=2) #adjusts acceptance probability rate


#### STEP 5: OPTIMIZE A SET OF MULTIPLEX PRIMERS BY MINIMIZING DIMER LOAD
# To run optimization once (simple problems only, multiple runs still recommended):
mw.optimizeMultiplex(PRIMER_FASTA = os.path.join(OUTDIR1, 'FilteredPrimers.fa'), 
                     DIMER_SUMS = DIMER_TOTS,
                     DIMER_TABLE = os.path.join(OUTDIR2, 'PrimerPairInteractions_wide.csv'), 
                     OUTPATH = os.path.join(OUTDIR3,"OUTNAME"), 
                     N_LOCI = N_LOCI, 
                     KEEPLIST = KEEPLIST_FA, #KEEPLIST_FA,
                     deltaG = DELTAG, #True for deltaG optimization, False for standard optimization
                     VERBOSE = False,#set to true to print dimers at each change
                     SIMPLE = 5000, # iterations for simple iterative improvement optimization (default=5000)
                     ITERATIONS = 1000, # iterations per simulated annealing cycle (default=1000) 
                     CYCLES = 10, # number of simualted annealing iterations to run (default=10)
                     BURNIN = 100, # iterations for sampling dimer cost space to adaptively set SA temps (default=100)
                     DECAY_RATE = 0.95, # temperature decay parameter for SA temps (default=0.98)
                          # closer to 1 - least conservative, explores more cost space at higher risk
                          # closer to 0 - most conservative, explores less cost space at lower risk
                          # recommendations: 0.90-0.98, higher with fewer iterations
                     T_INIT = None, # starting temp for fixed SA schedule (default=0.1)
                     T_FINAL = 0.01, # ending temp for fixed SA schedule (default=None, i.e., adaptively set based on costs observed in BURNIN)
                          # temperatures=0 is equivalent to simple iterative improvement, while 
                          # higher temperatures explore more of the cost space at higher risk of accepting dimers
                          # recommended initial fixed schedule is T_INIT~2 and T_FINAL=0.1
                     PROB_ADJ = 2,# adjusts dimer acceptance probabilities (default=2)
                          # increase if too many dimers are being accepted during simulated annealing, 
                          # decrease if local optima are not being overcome
                     SEED = None,#primer set from previous optimization run to start with, in CSV format
                     MAKEPLOT = False)#whether to run plotSAtemps within
# Outputs are found under {OUTDIR}/3_OptimizedMultiplexes/ and include:
# {OUTPATH}_primers.csv: primers included in the optimized multiplex
# OUTPATH_dimers.csv : pairwise dimer loads of primer pairs within the optimized multiplex
# OUTPATH_costsTrace.csv : Recorded change in dimer load at each accepted swap in the multiplex.
# OUTPATH_DimerLoad.png : Trace plot of cost function through algorithm progression

# Convert output from single run to FASTA
mw.CSVtoFASTA(IN_CSV = os.path.join(OUTDIR3, "OUTNAME_primers.csv"), 
              OUT_FA = os.path.join(OUTDIR3, "OUTNAME_primers.fasta"))


## The above function just runs the optimization process once. I highly recommend
## multiplex runs and then using the best option, especially for complex problems
## where outputs will vary. 
## The multipleOptimizations function will run optimizeMultiplex N_RUNS times:
mw.multipleOptimizations(N_RUNS = 10, #number of optimization runs
                         PRIMER_FA = os.path.join(OUTDIR1, 'FilteredPrimers.fa'), 
                         DIMER_SUMS = os.path.join(OUTDIR2, 'PrimerPairInteractions_sum.csv'), 
                         DIMER_TABLE = os.path.join(OUTDIR2, 'PrimerPairInteractions_wide.csv'), 
                         OUTPATH = os.path.join(OUTDIR3,"OUTNAME"), 
                         N_LOCI = N_LOCI,
                         deltaG = DELTAG,#False for standard optimization, True for deltaG optimization
                         KEEPLIST = KEEPLIST_FA, 
                         TIMEOUT = 10,#time allowed per run- runs 10 minutes will break
                         VERBOSE=False,#set to true to print dimers at each change
                         SIMPLE=5000, # iterations for simple iterative improvement optimization (default=5000)
                         ITERATIONS=1000, # iterations per simulated annealing cycle (default=1000) 
                         CYCLES=10, # number of simulated annealing cycles to run (default=10)
                         BURNIN=100, # iterations for sampling dimer cost space to adaptively set SA temps (default=100)
                         DECAY_RATE=0.95, # temperature decay parameter for SA temps (default=0.98)
                          # closer to 1 - least conservative, explores more cost space at higher risk
                          # closer to 0 - most conservative, explores less cost space at lower risk
                          # recommendations: 0.90-0.98, higher with fewer iterations
                         T_INIT=None, # starting temp for fixed SA schedule- higher: more hill-climbing
                         T_FINAL=0.01, # ending temp for fixed SA schedule (default=None, i.e., adaptively set based on costs observed in BURNIN)
                         PROB_ADJ=2,# adjusts dimer acceptance probabilities (default=2)
                          # increase if too many dimers are being accepted during simulated annealing, 
                          # decrease if local optima are not being overcome
                         SEED=None)#primer set from previous optimization run to start with, in CSV format
# Outputs will be saved to {OUTDIR}/3_OptimizedMultiplexes, including:
# - {OUTNAME}_RunSummary.csv : summarizes 
# - Run{x}_{OUTNAME}_primers.fa : Multiplexes with minimum dimer load across all runs
# - Final_Dimers: *_dimers.csv files for each run
# - Final_Primers: *_primers.csv files for each run
# - Plots_Dimer_Load: *_DimerLoad.png files for each run
# - Trace_Dimer_Load: *costsTrace.csv files for each run

## If you want to assess the output panel:
mw.assessPanel(os.path.join(OUTDIR3, "Run01_OUTNAME.fa"))


### IMPORTANT! ADDITIONAL SPECIFICITY CHECKS SHOULD BE RUN ON THIS PANEL
### (e.g., PRIMER-BLAST) BEFORE LAB TESTING!
