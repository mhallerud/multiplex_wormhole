#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: MULTIPLEX PRIMER DESIGN with multiplex_wormhole

Created on Tue Dec 19 19:50:45 2023
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
        
    WHITELIST_FA : A FASTA formatted file containing adapter-ligated primer sequences from a current multiplex
        assay that you are looking to add to. 
        Forward primers must be designated with ".FWD" or ".FW" as a suffix and reverse primers with ".REV"
    
    Note that locus IDs in the TEMPLATES file and primer IDs in the WHITELIST_FA must be unique
    for multiplex_wormhole to function properly. IDs may not contain periods.
"""


# load dependencies and modules
import os

# set working dir to location of multiplex_wormhole
os.chdir('/Users/maggiehallerud/Desktop/multiplex_wormhole')
from scripts.primer3_batch_design import main as primer3BatchDesign
from scripts.filter_primers import main as filterPrimers
from scripts.check_primer_specificity import main as specificityCheck
from scripts.add_whitelist_to_fasta import main as addWhitelistFasta
from scripts.tabulate_MFEprimer_dimers import main as tabulateDimers
from scripts.optimize_primers import main as optimizeMultiplex
from scripts.plot_SA_temps import main as plotSAtemps


## SET INPUTS:
os.chdir("/Users/maggiehallerud/Desktop/GrayFoxSNPs/insilico_design")
MFEprimer_PATH='/Users/maggiehallerud/Marten_Primer_Design/Plate1_First55Pairs_Sep2023/mfeprimer-3.2.7-darwin-10.6-amd64'
PRIMER3_PATH='/Users/maggiehallerud/primer3/src/primer3_core' #path to primer3 location
TEMPLATES='../Input_SNPs/GrayFox_primer3Templates.csv'
WHITELIST_FA=None#'/Users/maggiehallerud/Desktop/Marten_Fisher_Population_Genomics_Results/Marten/SNPpanel/Panel2_200initialpairs/Marten_Panel1_whitelist.fa'
OUTDIR='individualID'
#GENOME='/Users/maggiehallerud/Marten_Primer_Design/Plate2_Oct2023/0_Inputs/CoastalMartens.maf30.CENSORmask.fa'
N_LOCI = 75



## Step 0: Set up output directory structure & copy inputs to it
if not os.path.exists(OUTDIR):
    os.mkdir(OUTDIR)
INPUTDIR = os.path.join(OUTDIR, '0_Inputs')
if not os.path.exists(INPUTDIR):
    os.mkdir(INPUTDIR)
os.system("cp "+TEMPLATES+" "+INPUTDIR)
if WHITELIST_FA is not None:
    os.system("cp "+WHITELIST_FA+" "+INPUTDIR)
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
primer3BatchDesign(TEMPLATES, OUTDIR, PRIMER3_PATH)
# Outputs are found in the 1_InitialPrimers folder. there is a *.out and *.err file per locus



## Step 2: filter out primers with dimers
filterPrimers(PRIMER_DIR = os.path.join(OUTDIR, '1_InitialPrimers'), 
              OUTPATH = os.path.join(OUTDIR2,'FilteredPrimers'),
              Tm_LIMIT=45, 
              dG_HAIRPINS=-2000, 
              dG_END_LIMIT=-5000,
              dG_MID_LIMIT=-10000)
# Outputs are found under 2_FilteredPrimers/FilteredPrimers*



## Step 3A Checking primer specificity 
# Step 3A: Check primer specificity against provided loci
specificity_output = os.path.join(OUTDIR2,'SpecificityCheckTemplates')
specificityCheck(PRIMERS = os.path.join(OUTDIR2,'FilteredPrimers.csv'),
                 TARGET = TEMPLATES, 
                 OUTPATH = specificity_output)
# Outputs are found under 2_FilteredPrimers/SpecficityCheckTemplates*

INPUT=os.path.join(OUTDIR2, 'SpecificityCheckTemplates_passed.fa')#this is input for step 4



## IMPORTANT NOTE: If you have previous loci that you want included in this panel, now is the time to add them.
## BEFORE RUNNING THIS STEP: Check that whitelist IDs must have suffixes of ".FW" and ".REV", and may not contain any other periods.

## Here's a helper script to automate combining these files:
if WHITELIST_FA is not None:
    WHITELIST_FA = os.path.join(INPUTDIR, os.path.basename(WHITELIST_FA))
    if os.path.exists(WHITELIST_FA):
        addWhitelistFasta(specificity_output+'_passed.fa', WHITELIST_FA)
        # adjust input for step 4 to include whitelist
        INPUT = os.path.join(OUTDIR2, 'SpecificityCheckTemplates_passed_plusWhitelist.fa')

## Here is another helper script to convert CSV format primers to FA format:
#csvToFasta(IN_CSV, ID_FIELD, SEQ_FIELD, OUT_FA)




## Step 4: Predict primer dimers using MFEprimer
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
os.system(MFEprimer_PATH+" dimer -i "+INPUT+" -o "+END_DIMERS+" -d -5 -s 3 -m 70 --diva 3.8 --mono 50 --dntp 0.25 --oligo 50 -p")



## Step 5: Convert MFEprimer dimer report to table formats
## NOTE: This is the most computationally intensive step. 
## It will run substantially faster if you leave the 4th argument blank 
## (which means pairwise interactions between individual primers won't be calculated)
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
plotSAtemps(OUTPATH=os.path.join(OUTDIR4, 'TestingDefaults_75loci'),
            PRIMER_FASTA=os.path.join(OUTDIR2, 'SpecificityCheckTemplates_passed.fa'), 
            DIMER_SUMS=os.path.join(OUTDIR3, 'PrimerPairInteractions_sum.csv'), 
            DIMER_TABLE=os.path.join(OUTDIR3, 'PrimerPairInteractions_wide.csv'), 
            N_LOCI=75, #number of target loci in panel
            WHITELIST=WHITELIST_FA, 
            SEED=None, 
            BURNIN=100)#number iterations with dimer loads used to sample cost space
# decay rate closer to 1: 
plotSAtemps(OUTPATH=os.path.join(OUTDIR4, 'TestingSAparams_75loci_decayRate95'),
            # dimer counts to plot and calculate temps from (if not set)
            MIN_DIMER=1,
            MAX_DIMER=5, #update this value based on the max observed in the default plot!
            # parameter determining temperature decay in negative exponential
            DECAY_RATE=0.95, #default is 0.98
            # initial temperature to start from - higher=more risk accepted
            T_INIT=2, 
            # final temperature to stop at - 1=no risk, higher=more risk accepted
            T_FINAL=0, 
            #proportion of max dimer load considered when setting temperature schedule
            #closer to 0 = accepts fewer errors
            DIMER_ADJ=0.1,
            # adjustment for dimer acceptance probabilities- 1=no adjustment, higher values=lower dimer acceptance
            PROB_ADJ=1)



## Step 6: Design a set of multiplex primers by minimizing predicted dimer formation
# N_LOCI here is the number of loci you want in the final panel (including whitelist loci)
# To run once:
optimizeMultiplex(PRIMER_FASTA = os.path.join(OUTDIR2, 'SpecificityCheckTemplates_passed.fa'), 
                  DIMER_SUMS = os.path.join(OUTDIR3, 'PrimerPairInteractions_binary_sum.csv'), 
                  DIMER_TABLE = os.path.join(OUTDIR3, 'PrimerPairInteractions_binary_wide.csv'), 
                  OUTPATH = os.path.join(OUTDIR4,"Run01_75Loci_MAF30"), 
                  N_LOCI = N_LOCI, 
                  WHITELIST = WHITELIST_FA,
                  VERBOSE=True)
# Outputs are found under 4_OptimizedSets/*

## NOTE: I recommend rerunning this multiple times and taking the best option, since this is a 
## random process and each run may be slightly different.
## Here's a helper function for running N times:
from scripts.multiple_run_optimization import multipleOptimizations
multipleOptimizations(N_RUNS = 15, 
                      PRIMER_FA = os.path.join(OUTDIR2, 'SpecificityCheckTemplates_passed.fa'), 
                      DIMER_SUMS = os.path.join(OUTDIR3, 'PrimerPairInteractions_binary_sum.csv'), 
                      DIMER_TABLE = os.path.join(OUTDIR3, 'PrimerPairInteractions_binary_wide.csv'), 
                      OUTPATH = os.path.join(OUTDIR4,"MAF30_75loci"), 
                      N_LOCI = N_LOCI, 
                      WHITELIST = WHITELIST_FA, 
                      TIMEOUT = 120)#time allowed per run- runs >10 minutes will break

