#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: MULTIPLEX PRIMER DESIGN & ASSESSMENT with multiplex_wormhole
Purpose: multiplex_wormhole optimizes primer design for multiplex amplicon sequencing 
    by minimizing predicted pairwise dimers, and assesses existing multiplexes. 
    The target audience is for SNP panel development, however the process is transferable 
    to other targeted sequencing applications.
Dependencies: primer3-py (developed w/ v2.0.0)
              pandas (developed w/ v1.4.4)
              matplotlib (developed w/ v3.5.2)
              MFEprimer v3 --> install and setup with setup_mfeprimer.py
              
Input preparation:
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
    
Created on Tue Dec 19 19:50:45 2023
@author: maggiehallerud
"""

#### LOAD MULTIPLEX WORMHOLE FUNCTIONS & DEPENDENCIES ####
import sys
import glob


## CHANGE TO YOUR MULTIPLEX_WORMHOLE PATH!
## NO SPACES IN FILEPATHS ALLOWED HERE!
MULTIPLEX_WORMHOLE = "/Users/maggiehallerud/Desktop/multiplex_wormhole/"
sys.path.append(MULTIPLEX_WORMHOLE+"/src")

## INSTALL MFEprimer - THIS ONLY NEEDS TO BE RUN ONCE!
from setup_mfeprimer import main as setup_mfeprimer
setup_mfeprimer()
## IF THIS FAILS, YOU CAN ALSO MANUALLY DOWNLOAD & SET THE PATH TO MFEPRIMER BELOW:
## (Remember- no spaces in filepaths)
## If downloading manually, you may also need to change permission. In terminal, run e.g.:
## chmod +x /yourpath/multiplex_wormhole/src/mfeprimer*
MFEprimer_PATH=glob.glob(MULTIPLEX_WORMHOLE+"src/*mfeprimer*")[0]
#MFEprimer_PATH='/Users/maggiehallerud/Desktop/multiplex_wormhole/src/mfeprimer-3.2.7-darwin-10.6-amd64'


# load multiplex wormhole functions
import importlib
import os
import shutil
import subprocess
from scripts.primer3_batch_design import main as primer3BatchDesign
from scripts.tabulate_MFEprimer_dimers import main as tabulateDimers
from scripts.optimize_multiplex import main as optimizeMultiplex
from scripts.multiple_run_optimization import main as multipleOptimizations
from scripts.CSVtoFasta import main as CSVtoFASTA
plotASAtemps = importlib.import_module("plot_ASA_temps")




#### ALTERNATIVE WORKFLOW: PANEL ASSESSMENT ####
from panel_assessment import main as assessPanel
# INPUT: FASTA or CSV (PrimerID, Sequence) file with primers, PrimerIDs following rules
# described above (e.g., "MACA1.FWD" & "MACA1.REV")
assessPanel("Primers.fasta")
# OUTPUTS: # primer pairs, total # pairwise dimers, # primer pairs forming dimers
# also dimer files from MFEprimer and dimer tables from tabulateDimers

    
    
#### PRIMARY WORKFLOW: OPTIMIZED PANEL DESIGN FOR MULTIPLEX PCR ####
## SET INPUTS:
os.chdir("MW_TEST/")#path to project folder
TEMPLATES = "/multiplex_wormhole/examples/Input_Templates.csv"#CSV containing candidate sequences (path relative to project folder)
KEEPLIST_FA = None #"MartenPanel1.fa" #FASTA containing previously designed primer set
OUTDIR = "TEST" # folder name where outputs will be saved
N_LOCI = 50 # target panel size (# sequences amplified)
DELTAG = False #set to True if you want to use deltaG optimization algorithm


## Step 0: Set up output directory structure & copy inputs to it
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

    
## Step 1: batch design of primers
# NOTE: This script includes all of the primer design settings!
# These can be adjusted directly in the primer3_batch_design.py script
# or by providing a dictionary to SETTINGS.
# details on settings: https://htmlpreview.github.io/?https://github.com/primer3-org/primer3/blob/v2.6.1/src/primer3_manual.htm#globalTags 
primer3BatchDesign(TEMPLATES, 
                   os.path.join(OUTDIR1, "FilteredPrimers"), 
                   Tm_LIMIT = 45,
                   dG_HAIRPINS = -2,
                   dG_END_LIMIT = -4, # lower 
                   dG_MID_LIMIT = -8, # lower limit for deltaG of all other dimers
                   KEEPLIST = KEEPLIST_FA, #FASTA for keeplist
                   ENABLE_BROAD = False, #use broader settings if no primers for template?
                   SETTINGS = None) #primer3 settings in dictionary {} format 
# Outputs are 1_PrimerDesign folder and include a FASTA & CSV with primer details.


## Step 2: Predict primer dimers using MFEprimer
# Set input based on whether keeplist is provided or not
if KEEPLIST_FA is None:
    INPUT = os.path.join(OUTDIR1, "FilteredPrimers.fa")
else:
    INPUT = os.path.join(OUTDIR1, "FilteredPrimers_plusKeeplist.fa")

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
subprocess.call(MFEprimer_PATH+" dimer -i "+INPUT+" -o "+ALL_DIMERS+" -d -8 -s 3 -m 50 --diva 3.8 "+
          "--mono 50 --dntp 0.25 --oligo 50", shell=True)
subprocess.call(MFEprimer_PATH+" dimer -i "+INPUT+" -o "+END_DIMERS+" -d -4 -s 3 -m 70 --diva 3.8 "+
          "--mono 50 --dntp 0.25 --oligo 50 -p", shell=True)



## Step 3: Convert MFEprimer dimer report to table formats
## NOTE: This is the most computationally intensive step. 
## It will run substantially faster if you leave the 4th argument blank 
## (which means pairwise interactions between individual primers won't be calculated)
tabulateDimers(ALL_DIMERS, 
               END_DIMERS, 
               os.path.join(OUTDIR2, 'PrimerPairInteractions'), 
               "False",#os.path.join(OUTDIR3, 'RawPrimerInteractions'))#specify this parameter if you care about per-primer dimers (Rather than just sums per primer pair)
               DELTAG)
# Outputs are found under 2_PredictedDimers/PrimerPairInteractions*
# names will depend on deltaG....
if DELTAG:
    DIMER_TOTS = os.path.join(OUTDIR2, 'PrimerPairInteractions_mean.csv')
else:
    DIMER_TOTS = os.path.join(OUTDIR2, 'PrimerPairInteractions_sum.csv')



## Step 4: Explore temperature space for simulated annealing
## There are two ways to run this script: one calculates temperatures and dimer loads based on the problem at hand, 
## the other uses pre-specified temperatures and dimer loads.
## I recommend first running using files from the problem, then using the values observed in the outputs to explore 
## parameters around the defaults.
plotASAtemps.main(OUTPATH=os.path.join(OUTDIR3, 'TestingASAparams_defaults'),
                  PRIMER_FASTA=os.path.join(OUTDIR1, 'FilteredPrimers.fa'), 
                  DIMER_SUMS=DIMER_TOTS,
                  DIMER_TABLE=os.path.join(OUTDIR2, 'PrimerPairInteractions_wide.csv'), 
                  N_LOCI=N_LOCI, #number of target loci in panel
                  KEEPLIST=KEEPLIST_FA, 
                  SEED=None, #this would be an output from optimizeMultiplex
                  BURNIN=100,#number iterations with dimer loads used to sample cost space
                  deltaG=DELTAG)
# decay rate closer to 1: 
plotASAtemps.main(OUTPATH=os.path.join(OUTDIR3, 'TestingASAparams_decayRate98'),
                  # dimer counts to plot and calculate temps from (if not set)
                  MIN_DIMER=1,
                  MAX_DIMER=5, #update this value based on the max observed in the default plot!
                  # parameter determining temperature decay in negative exponential
                  DECAY_RATE=0.98, #default is 0.95
                  # initial temperature to start from - higher=more risk accepted
                  T_INIT=2, 
                  # final temperature to stop at - 1=no risk, higher=more risk accepted
                  T_FINAL=0, #DEFAULT: 0.1
                  #proportion of max dimer load considered when setting temperature schedule
                  #closer to 0 = accepts fewer errors
                  DIMER_ADJ=0.1,
                  # adjustment for dimer acceptance probabilities- 1=no adjustment, higher values=lower dimer acceptance
                  PROB_ADJ=1)#DEFAULT=2


## Step 5: Design a set of multiplex primers by minimizing predicted dimer formation
# N_LOCI here is the number of loci you want in the final panel (including keeplist loci)
# To run once:
optimizeMultiplex(PRIMER_FASTA = os.path.join(OUTDIR1, 'FilteredPrimers.fa'), 
                  DIMER_SUMS = DIMER_TOTS,
                  DIMER_TABLE = os.path.join(OUTDIR2, 'PrimerPairInteractions_wide.csv'), 
                  OUTPATH = os.path.join(OUTDIR3,"OUTNAME"), 
                  N_LOCI = N_LOCI, 
                  KEEPLIST = None, #KEEPLIST_FA,
                  deltaG = DELTAG, #True for deltaG optimization, False for standard optimization
                  VERBOSE=False,#set to true to print dimers at each change
                  SIMPLE=3000, # iterations for simple iterative improvement optimization (default=5000)
                  ITERATIONS=5000, # iterations for simulated annealing optimization (default=10000) 
                  BURNIN=100, # iterations for sampling dimer cost space to adaptively set SA temps (default=100)
                  DECAY_RATE=0.98, # temperature decay parameter for SA temps (default=0.98)
                      # closer to 1 - least conservative, explores more cost space at higher risk
                      # closer to 0 - most conservative, explores less cost space at lower risk
                      # recommendations: 0.90-0.98, higher with fewer iterations
                  T_INIT=None, # starting temp for fixed SA schedule (default=0.1)
                  T_FINAL=None, # ending temp for fixed SA schedule (default=None, i.e., adaptively set based on costs observed in BURNIN)
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
                  SEED=None,#primer set from previous optimization run to start with, in CSV format
                  MAKEPLOT=False)#whether to run plotSAtemps within
# Outputs are found under 3_OptimizedSets/* and include:
# OUTPATH_primers.fa & _primers.csv: primers included in the optimized multiplex
# OUTPATH_dimers.fa : dimer load of each primer pair within the optimized multiplex
# OUTPATH_ASA_costs.csv : Recorded change in dimer load at each accepted swap in the multiplex.
# OUTPATH_dimerload.png : Trace plot of cost function through algorithm progression


## The above function just runs the optimization process once. I highly recommend
## multiplex runs and then using the best option, especially for complex problems
## where outputs will vary. 
## The multipleOptimizations function will run optimizeMultiplex N_RUNS times:
multipleOptimizations(N_RUNS = 10, 
                      PRIMER_FA = os.path.join(OUTDIR1, 'FilteredPrimers.fa'), 
                      DIMER_SUMS = os.path.join(OUTDIR2, 'PrimerPairInteractions_sum.csv'), 
                      DIMER_TABLE = os.path.join(OUTDIR2, 'PrimerPairInteractions_wide.csv'), 
                      OUTPATH = os.path.join(OUTDIR3,"OUTNAME"), 
                      N_LOCI = 100, 
                      deltaG = DELTAG,#False for standard optimization, True for deltaG optimization
                      KEEPLIST = KEEPLIST_FA, 
                      TIMEOUT = 10,#time allowed per run- runs 10 minutes will break
                      VERBOSE=False,#set to true to print dimers at each change
                      SIMPLE=3000, # iterations for simple iterative improvement optimization (default=5000)
                      ITERATIONS=5000, # iterations for simulated annealing optimization (default=10000) 
                      BURNIN=100, # iterations for sampling dimer cost space to adaptively set SA temps (default=100)
                      DECAY_RATE=0.98, # temperature decay parameter for SA temps (default=0.98)
                          # closer to 1 - least conservative, explores more cost space at higher risk
                          # closer to 0 - most conservative, explores less cost space at lower risk
                          # recommendations: 0.90-0.98, higher with fewer iterations
                      T_INIT=2, # starting temp for fixed SA schedule (default=0.1)
                      T_FINAL=0.1, # ending temp for fixed SA schedule (default=None, i.e., adaptively set based on costs observed in BURNIN)
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
#OUTPUTS: 3_OptimizedSets/{OUTNAME}_RunSummary.csv = run summary
#3_OptimizedSets/{OUTNAME}_RunSummary.csv = run summary


## STEP 6: Convert selected primer set to FASTA format for additional screening
# CHANGE THESE TO THE "BEST" Run based on RunSummary.csv output!
CSVtoFASTA(IN_CSV = os.path.join(OUTDIR3, "OUTNAME_50loci_Run1_primers.csv"), 
           OUT_FA = os.path.join(OUTDIR3, "OUTNAME_50loci_Run1_primers.fasta"),
           ID_FIELD = "PrimerID",  
           SEQ_FIELD = "Sequence")
