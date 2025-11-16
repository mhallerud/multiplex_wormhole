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
import os
import sys
import datetime
import pandas as pd

# load multiplex wormhole functions
sys.path.append(os.path.dirname(__file__))
from scripts.primer3_batch_design import main as primer3BatchDesign
from scripts.filter_primers import main as filterPrimers
from scripts.check_primer_specificity import main as specificityCheck
from scripts.add_keeplist_to_fasta import main as addKeeplistFasta
from scripts.tabulate_MFEprimer_dimers import main as tabulateDimers
#from scripts.optimize_primers import main as optimizeMultiplex
#from scripts.plot_SA_temps import main as plotSAtemps
from scripts.multiple_run_optimization import multipleOptimizations
from scripts.extras.CSVtoFasta import main as CSVtoFASTA




def main(TEMPLATES, N_LOCI, OUTDIR, PREFIX=None, KEEPLIST_FA=None, N_RUNS=10, ITERATIONS=5000, SIMPLE=2000):
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
    Returns
    -------
    1. Designs primers for each candidate sequences in TEMPLATES
    2. Filtered primer sets in CSV format
    3. Pairwise dimers predicted between all filtered primer pairs
    4. Optimized multiplex primer set for target N_LOCI
    5. Predicted dimer loads for primers included in optimized multiplex
    6. Plots of simulated annealing temperature schedule and trace of dimer load
    """
    ## SET PATHS TO DEPENDENCIES
    MFEprimer_PATH='/Users/maggiehallerud/Marten_Primer_Design/Plate1_First55Pairs_Sep2023/mfeprimer-3.2.7-darwin-10.6-amd64'#full path to mfeprimer location
    PRIMER3_PATH='/Users/maggiehallerud/primer3/src/primer3_core' #full path to primer3 location

    
    ## Step 0: Set up output directory structure & copy inputs to it
    print("-----SETTING UP OUTPUT DIRECTORY STRUCTURE------")
    if not os.path.exists(OUTDIR):
        os.mkdir(OUTDIR)
    INPUTDIR = os.path.join(OUTDIR, '0_Inputs')
    if not os.path.exists(INPUTDIR):
        os.mkdir(INPUTDIR)
    os.system("cp "+TEMPLATES+" "+INPUTDIR)
    if KEEPLIST_FA is not None:
        os.system("cp "+KEEPLIST_FA+" "+INPUTDIR)
    OUTDIR2 = os.path.join(OUTDIR, "2_FilteredPrimers")
    if not os.path.exists(OUTDIR2):
        os.mkdir(OUTDIR2)
    OUTDIR3 = os.path.join(OUTDIR, "3_PredictedDimers")
    if not os.path.exists(OUTDIR3):
        os.mkdir(OUTDIR3)
    OUTDIR4 = os.path.join(OUTDIR, "4_OptimizedSets")
    if not os.path.exists(OUTDIR4):
        os.mkdir(OUTDIR4)
    
    # set suffix to current datetime if not given
    if PREFIX is None:
        PREFIX = str(datetime.datetime.now()).replace(" ","_").replace(".","_")
        
    ## Step 1: batch design of primers
    print("")
    print("-----BATCH DESIGNING PRIMERS------")
    primer3BatchDesign(TEMPLATES, OUTDIR, PRIMER3_PATH)
    # Outputs are found in the 1_InitialPrimers folder. there is a *.out and *.err file per locus
    
    
    
    ## Step 2: filter out primers with dimers
    print("")
    print("-----REMOVING PRIMERS WITH INTRA-PAIR DIMERS-----")
    filterPrimers(PRIMER_DIR = os.path.join(OUTDIR, '1_InitialPrimers'), 
                  OUTPATH = os.path.join(OUTDIR2,'FilteredPrimers'),
                  Tm_LIMIT=45, 
                  dG_HAIRPINS=-2000, 
                  dG_END_LIMIT=-5000,
                  dG_MID_LIMIT=-10000)
    # Outputs are found under 2_FilteredPrimers/FilteredPrimers*
    
    
    
    ## Step 3: Checking primer specificity 
    print("")
    print("-----CHECKING PRIMER SPECIFICITY AGAINST OTHER TEMPLATES-----")
    specificity_output = os.path.join(OUTDIR2,'SpecificityCheckTemplates')
    specificityCheck(PRIMERS = os.path.join(OUTDIR2,'FilteredPrimers.csv'),
                     TARGET = TEMPLATES, 
                     OUTPATH = specificity_output)
    # Outputs are found under 2_FilteredPrimers/SpecficityCheckTemplates*
    
    INPUT=specificity_output+"_passed.fa" #this is input for step 4
    
    
    ## IMPORTANT NOTE: If you have previous loci that you want included in this panel, now is the time to add them.
    ## BEFORE RUNNING THIS STEP: Check that keeplist IDs must have suffixes of ".FW" and ".REV", and may not contain any other periods.
    
    ## Here's a helper script to automate combining these files:
    if KEEPLIST_FA is not None:
        print("")
        print("-----ADDING KEEPLIST TO INPUT FASTA-----")
        KEEPLIST_FA = os.path.join(INPUTDIR, os.path.basename(KEEPLIST_FA))
        if os.path.exists(KEEPLIST_FA):
            addKeeplistFasta(specificity_output+'_passed.fa', KEEPLIST_FA)
            # adjust input for step 4 to include keeplist
            INPUT = os.path.join(OUTDIR2, 'SpecificityCheckTemplates_passed_plusKeeplist.fa')
    
    ## Here is another helper script to convert CSV format primers to FA format:
    #from scripts.extras.CSVtoFasta import main as CSVtoFASTA
    #csvToFasta(IN_CSV, ID_FIELD, SEQ_FIELD, OUT_FA)
    
    
    
    ## Step 4: Predict primer dimers using MFEprimer
    print("")
    print("-----RUNNING DIMER PREDICTION WITH MFEPRIMER-----")
   # NOTE: Originally, primers were checked via the PrimerSuite PrimerDimer function (http://www.primer-dimer.com/)
    # PrimerSuite PrimerDimerReport files can be converted to the necessary table/sum files using scripts/translate_primerSuite_report.R
    # I decided to transition to MFEprimer because primer-dimer.com returned an unreasonable number of dimers
    # set output paths
    ALL_DIMERS=os.path.join(OUTDIR3, 'MFEprimerDimers.txt')
    END_DIMERS=os.path.join(OUTDIR3, 'MFEprimerDimers_ends.txt')
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
    os.system(MFEprimer_PATH+" dimer -i "+INPUT+" -o "+ALL_DIMERS+" -d -8 -s 3 -m 50 --diva 3.8 --mono 50 --dntp 0.25 --oligo 50")
    os.system(MFEprimer_PATH+" dimer -i "+INPUT+" -o "+END_DIMERS+" -d -4 -s 3 -m 70 --diva 3.8 --mono 50 --dntp 0.25 --oligo 50 -p")
    
    
    
    ## Step 5: Convert MFEprimer dimer report to table formats
    ## NOTE: This is the most computationally intensive step. 
    ## It will run substantially faster if you leave the 4th argument blank 
    ## (which means pairwise interactions between individual primers won't be calculated)
    print("")
    print("-----TABULATING PREDICTED DIMERS-----")
    tabulateDimers(ALL_DIMERS, 
                   END_DIMERS, 
                   os.path.join(OUTDIR3, 'PrimerPairInteractions'), 
                   "False")#os.path.join(OUTDIR3, 'RawPrimerInteractions'))#specify this parameter if you care about per-primer dimers (Rather than just sums per primer pair)
    # Outputs are found under 3_PredictedDimers/PrimerPairInteractions*
    
    
    
    ## Step 6A: Explore temperature space for simulated annealing
    ## There are two ways to run this script: one calculates temperatures and dimer loads based on the problem at hand, 
    ## the other uses pre-specified temperatures and dimer loads.
    ## I recommend first running using files from the problem, then using the values observed in the outputs to explore 
    ## parameters around the defaults.
    # plotSAtemps(OUTPATH=os.path.join(OUTDIR4, 'TestingDefaults_50loci'),
    #             PRIMER_FASTA=os.path.join(OUTDIR2, 'SpecificityCheckTemplates_passed.fa'), 
    #             DIMER_SUMS=os.path.join(OUTDIR3, 'PrimerPairInteractions_sum.csv'), 
    #             DIMER_TABLE=os.path.join(OUTDIR3, 'PrimerPairInteractions_wide.csv'), 
    #             N_LOCI=N_LOCI, #number of target loci in panel
    #             KEEPLIST=KEEPLIST_FA, 
    #             SEED=None, #this would be an output from optimizeMultiplex
    #             BURNIN=100)#number iterations with dimer loads used to sample cost space
    # # decay rate closer to 1: 
    # plotSAtemps(OUTPATH=os.path.join(OUTDIR4, 'TestingSAparams_75loci_decayRate95'),
    #             # dimer counts to plot and calculate temps from (if not set)
    #             MIN_DIMER=1,
    #             MAX_DIMER=5, #update this value based on the max observed in the default plot!
    #             # parameter determining temperature decay in negative exponential
    #             DECAY_RATE=0.95, #default is 0.98
    #             # initial temperature to start from - higher=more risk accepted
    #             T_INIT=2, 
    #             # final temperature to stop at - 1=no risk, higher=more risk accepted
    #             T_FINAL=0, 
    #             #proportion of max dimer load considered when setting temperature schedule
    #             #closer to 0 = accepts fewer errors
    #             DIMER_ADJ=0.1,
    #             # adjustment for dimer acceptance probabilities- 1=no adjustment, higher values=lower dimer acceptance
    #             PROB_ADJ=1)
    
    
    ## Step 6: Design a set of multiplex primers by minimizing predicted dimer formation
    # N_LOCI here is the number of loci you want in the final panel (including keeplist loci)
 
    ## NOTE: I recommend rerunning this multiple times and taking the best option, since this is a 
    ## random process and each run may be slightly different.
    print("")
    print("-----STARTING OPTIMIZATION FOR MULTIPLEX PRIMER SET-----")
    multipleOptimizations(N_RUNS = N_RUNS, 
                          PRIMER_FA = os.path.join(OUTDIR2, 'SpecificityCheckTemplates_passed.fa'), 
                          DIMER_SUMS = os.path.join(OUTDIR3, 'PrimerPairInteractions_binary_sum.csv'), 
                          DIMER_TABLE = os.path.join(OUTDIR3, 'PrimerPairInteractions_binary_wide.csv'), 
                          OUTPATH = os.path.join(OUTDIR4, PREFIX), 
                          N_LOCI = N_LOCI, 
                          KEEPLIST = KEEPLIST_FA, 
                          TIMEOUT = 360,#time allowed per run- runs 30 minutes will break
                          VERBOSE=False,#set to true to print dimers at each change
                          SIMPLE=SIMPLE, # iterations for simple iterative improvement optimization (default=5000)
                          ITERATIONS=ITERATIONS, # iterations for simulated annealing optimization (default=10000) 
                          BURNIN=100, # iterations for sampling dimer cost space to adaptively set SA temps (default=100)
                          DECAY_RATE=0.98, # temperature decay parameter for SA temps (default=0.98)
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
                          SEED=None,#primer set from previous optimization run to start with, in CSV format
                          MAKEPLOT=False)#whether to run plotSAtemps during optimization
    #OUTPUT: MAF30_150loci_RunSummary.csv
    
    
    ## STEP 7: Convert selected primer set to FASTA format for additional screening
    print("")
    print("-----CONVERTING BEST PRIMER SETS INTO FASTAs FOR ADDITIONAL SCREENING-----")
    runs = pd.read_csv(os.path.join(OUTDIR4,PREFIX+"_RunSummary.csv"))
    run = runs['Run']
    run = [str(run[x]).zfill(2) for x in range(len(run))]
    dimers = runs['TotalDimers']
    print("The BEST multiplex had "+str(min(dimers))+ "total predicted dimers.")
    for i in range(len(dimers)):
        print("Run "+str(run[i])+" had "+str(dimers[i])+" dimers.")
        if dimers[i]==min(dimers):
            print(".....Converting to FASTA for additional screening")
            CSVtoFASTA(IN_CSV = os.path.join(OUTDIR4,PREFIX+"_Run"+run+"_primers.csv"), 
                       OUT_FA = os.path.join(OUTDIR4,"Run"+run+"_"+PREFIX+"_primers.fasta"))



if __name__=="__main__":
    if len(sys.argv) != 3:
        print("Multiplex wormhole takes 6 arguments:")
        print("   multiplex_wormhole.py TEMPLATES N_LOCI OUTDIR KEEPLIST_FA N_RUNS ITERATIONS=5000 SIMPLE=2000")
        print("defaults:    multiplex_wormhole.py TEMPLATES N_LOCI OUTDIR 'None' 10 5000 2000")
    else:
        main(sys.argv[1],
             sys.argv[2],
             sys.arv[3],
             sys.argv[4],
             sys.argv[5],
             sys.argv[6])