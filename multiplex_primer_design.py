#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 19:50:45 2023

@author: maggiehallerud

Dependencies: 1) primer3
              2) MFEprimer must be installed 
"""


# load modules
import os
import sys
import signal
from scripts.primer3_batch_design import main as primer3_batch_design
from scripts.filter_primers import main as filter_primers
from scripts.check_primer_specificity import main as specificity_check
from scripts.tabulate_MFEprimer_dimers import main as tabulate_dimers
from scripts.optimize_primers import main as optimize_multiplex


# COMMAND-LINE INPUTS
MFEprimer_PATH=sys.argv[1]
OUTDIR=sys.arv[2]
TEMPLATES=sys.argv[3]
GENOME=sys.argv[4]
WHITELIST_FA=sys.argv[5]
N_LOCI = sys.argv[6]
N_ITERATIONS=sys.argv[7]
N_RUNS=sys.argv[8]

## INPUTS:
MFEprimer_PATH='/Users/maggiehallerud/Marten_Primer_Design/Plate1_First55Pairs_Sep2023/mfeprimer-3.2.7-darwin-10.6-amd64'
#PRIMER3_PATH#path to primer3 location
OUTDIR='/Users/maggiehallerud/Marten_Primer_Design/Plate1_First55Pairs_Sep2023/TEST'
TEMPLATES='/Users/maggiehallerud/Marten_Primer_Design/Plate1_First55Pairs_Sep2023/MartenTemplates_MAF30-repBaseFilter-random300.csv'
GENOME='/Users/maggiehallerud/Marten_Primer_Design/Plate2_Oct2023/0_Inputs/CoastalMartens.maf30.CENSORmask.fa'
WHITELIST_FA='/Users/maggiehallerud/Marten_Primer_Design/Plate1_First55Pairs_Sep2023/TEST/WHITELIST_FAKE.fa'
N_LOCI = 50
N_ITERATIONS=5000
N_RUNS=10



def main():
    ## Set up output directory structure
    if not os.path.exists(OUTDIR):
        os.mkdir(OUTDIR)
    OUTDIR2 = os.path.join(OUTDIR, "2_FilteredPrimers")
    if not os.path.exists(OUTDIR2):
        os.mkdir(OUTDIR2)
    OUTDIR3 = os.path.join(OUTDIR, "3_PredictedDimers")
    if not os.path.exists(OUTDIR3):
        os.mkdir(OUTDIR3)
    OUTDIR4 = os.path.join(OUTDIR, "4_OptimizedSets")
    if not os.path.exists(OUTDIR4):
        os.mkdir(OUTDIR4)
    
    
    ## Step 1: batch design of primers
    primer3_batch_design(TEMPLATES, OUTDIR)
    # for 200 loci: 45.5 seconds
    
    
    ## Step 2: filter out primers with dimers
    filter_primers(os.path.join(OUTDIR, '1_InitialPrimers'), 
                   os.path.join(OUTDIR2,'FilteredPrimers'), 
                   Tm_LIMIT=45, 
                   dG_HAIRPINS=-2000, 
                   dG_END_LIMIT=-5000,
                   dG_MID_LIMIT=-10000)
    # for 200 loci: 0.51 seconds
    
    
    ## Step 3A Checking primer specificity 
    # Step 3A: Check primer specificity against provided loci
    specificity_check(os.path.join(OUTDIR2,'FilteredPrimers.csv'), 
                      TEMPLATES, 
                      os.path.join(OUTDIR2,'SpecificityCheckTemplates'))
    # for 200 loci vs 200 loci: 6.15 seconds
    
    # Step 3B: Check primer specificity against a genome
    specificity_check(os.path.join(OUTDIR2, 'SpecificityCheckTemplates_passed.csv'), 
                      GENOME, 
                      os.path.join(OUTDIR2,'SpecificityCheckGenome'))
    # for 200 loci vs 3160 reads: 10.8 seconds
    
    
    
    ## IMPORTANT NOTE: If you have previous loci that you want included in this panel, now is the time to add them.
    ## This can be done by adding these loci into the SpecificityCheckGenome_passed.fa file.
    ## There is a helper script to convert CSV format primers to FA format:
    #csvToFasta(IN_CSV, ID_FIELD, SEQ_FIELD, OUT_FA)
    
    
    
    ## Step 4: Predict primer dimers using MFEprimer
    # NOTE: Originally, primers were checked via the PrimerSuite PrimerDimer function (http://www.primer-dimer.com/)
    # PrimerSuite PrimerDimerReport files can be converted to the necessary table/sum files using scripts/translate_primerSuite_report.R
    # I decided to transition to MFEprimer because primer-dimer.com returned an unreasonable number of dimers
    # set input / output paths
    INPUT=os.path.join(OUTDIR2, 'SpecificityCheckGenome_passed.fa')
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
    os.system(MFEprimer_PATH+" dimer -i "+INPUT+" -o "+ALL_DIMERS+" -d -6 -s 3 -m 40 --diva 3.8 --mono 50 --dntp 0.25 --oligo 50")
    os.system(MFEprimer_PATH+" dimer -i "+INPUT+" -o "+END_DIMERS+" -d -2.5 -s 3 -m 40 --diva 3.8 --mono 50 --dntp 0.25 --oligo 50 -p")
    # for 200 loci: 2.6 seconds (step 1: 1.6 sec, step 2: 0.9 sec)
    
    
    ## Step 5: Convert MFEprimer dimer report to table formats
    ## NOTE: This is the most computationally intensive step. 
    ## It will run substantially faster if you leave the 4th argument blank 
    ## (which means pairwise interactions between individual primers won't be calculated)
    tabulate_dimers(ALL_DIMERS, 
                    END_DIMERS, 
                    os.path.join(OUTDIR3, 'PrimerPairInteractions'), 
                    "False")#os.path.join(OUTDIR3, 'RawPrimerInteractions'))
    # for 200 loci (primer pairs and individual primers): 3 minutes 54 seconds
    # for 200 loci (primer pairs only): 48 seconds
    
    
    ## Step 6: Design a set of multiplex primers by minimizing predicted dimer formation
    ## NOTE: I recommend rerunning this multiple times and taking the best option, since this is a 
    ## random process and each run may be different.
    if WHITELIST_FA is None:
        for run in range(N_RUNS):
            signal.alarm(300) # set timer for 5 minutes- this helps to skip infinite loops in optimization process
            try:
                optimize_multiplex(PRIMER_FASTA=os.path.join(OUTDIR2, 'SpecificityCheckGenome_passed.fa'), 
                                   DIMER_SUMS=os.path.join(OUTDIR3, 'PrimerPairInteractions_binary_sum.csv'), 
                                   DIMER_TABLE=os.path.join(OUTDIR3, 'PrimerPairInteractions_binary_wide.csv'), 
                                   OUTPATH=os.path.join(OUTDIR4,"Run"+str(run+1)), 
                                   N_LOCI=N_LOCI, 
                                   ITERATIONS=N_ITERATIONS, 
                                   WHITELIST=None)
            # go to top of for loop if timed out
            except TimeoutException:
                continue
            # reset alarm
            else:
                signal.alarm(0)                
            print(" ")
    else:
        for run in range(N_RUNS):
            signal.alarm(300) # set timer for 5 minutes- this helps to skip infinite loops in optimization process
            try:
                optimize_multiplex(PRIMER_FASTA=os.path.join(OUTDIR2, 'SpecificityCheckGenome_passed.fa'), 
                                   DIMER_SUMS=os.path.join(OUTDIR3, 'PrimerPairInteractions_binary_sum.csv'),
                                   DIMER_TABLE=os.path.join(OUTDIR3, 'PrimerPairInteractions_binary_wide.csv'), 
                                   OUTPATH=os.path.join(OUTDIR4,"Run"+str(run+1)), 
                                   N_LOCI=N_LOCI, 
                                   ITERATIONS=N_ITERATIONS, 
                                   WHITELIST=WHITELIST_FA)
            except TimeoutException:
                continue
            #reset alarm
            else:
                signal.alarm(0)
            print(" ")


### Set up timeout exception behavior
class TimeoutException(Exception):   # Custom exception class
    pass

def timeout_handler(signum, frame):   # Custom signal handler
    raise TimeoutException

# Change the behavior of SIGALRM
signal.signal(signal.SIGALRM, timeout_handler)


if __name__=="__main__":
    main()
