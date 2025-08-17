#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: BATCH PRIMER DESIGN TOOL
Last edited on Tue Dec 19 18:28:09 2023
@author: maggiehallerud

Purpose: The batch primer design tool automates primer design for many loci. Primers are designed via
        primer3 based on default settings in primer3_BaseSettings.txt and backup settings primer3_BroadSettings.txt
Dependencies: primer3 must be installed and the filepath to primer3_core must be updated on line 19 in the primer3.sh file
Inputs: 1) a CSV file with the locus ID in column 1, DNA sequence in column 2, and target amplicon (format: start position,length)
        This file can be created using openPrimeR, see lines 1-90 in the helper script "0_create_in_templates.R"
        2) an output directory where primer details are saved
Output: primer3 file outputs containing information on primer pairs including sequences, annealing temperatures, and predicted 
        secondary structurees
"""



# load dependencies
import os
import sys
import shutil
import csv
import glob



def main(IN_CSV, OUTDIR, PRIMER3_PATH):
    """
    IN_CSV : csv filepath
        field1=ID, field2=DNA sequence, field3=amplicon (start_bp,length)

    OUTDIR : directory path
        an output directory where primer details are saved
    
    PRIMER3_PATH : primer3_core path
    	location of primer3_core software
    ------
    Designs primers for all loci in IN_CSV; 
    Returns primer output files in OUTDIR/1_InitialPrimers
    """
    # find script directory 
    SCRIPTPATH=os.path.dirname(__file__)

    # set up output folder structure
    if not os.path.exists(OUTDIR):
        os.mkdir(OUTDIR)
    
    # set up output directory to save primers
    outprimers = os.path.join(OUTDIR,'1_InitialPrimers')
    # if the folder already exists, delete and remake
    # this avoids issues building up between runs
    if os.path.exists(outprimers):
        shutil.rmtree(outprimers)
    os.mkdir(outprimers)
    
    # define paths to scripts and settings
    primer3_sh = os.path.join(SCRIPTPATH, 'primer3.sh')
    basedir = os.path.dirname(SCRIPTPATH)
    strict = os.path.join(basedir, 'primer3_settings/primer3_Base_NoSecondaryFilters.txt')
    broad = os.path.join(basedir, 'primer3_settings/primer3_Broad_NoSecondaryFilters.txt')
    
    # load templates CSV file with SEQUENCE_ID, SEQUENCE_TEMPLATE, SEQUENCE_TARGET in columns
    # use create_in_templates.R script lines 1-90 to create this CSV
    print("Reading in sequences......")
    templates = [] # empty array to store lines
    with open(IN_CSV, 'r', newline="\n") as file:
        reader = csv.reader(file, delimiter=",")
        next(reader) # skip header line
        for line in reader: 
            templates.append(line)
    print("     "+str(len(templates))+" target sequences input")
    
    # design FW and REV primers for each sequence
    print("Designing initial primers......")
    for row in range(len(templates)):
        # pull in inputs
        ids = templates[row][0]
        seq = templates[row][1]
        snp = templates[row][2]
        
        # run primer3 based on strict settings
        #print(" - Designing primers for "  + ids)
        os.system(primer3_sh +' '+ PRIMER3_PATH +' '+ ids +' '+ seq +' '+ snp +' '+ strict +' '+ outprimers)
    
    	# rerun with broad settings if no primers were found
        primers = glob.glob(os.path.join(outprimers, ids +'.out'))
        if len(primers)==0:
            print("     No primers found for "+ids+"! Retrying with broader settings.")
            os.system(primer3_sh +' '+ PRIMER3_PATH + ' '+ ids +' '+ seq +' '+ snp +' '+ broad +' '+ outprimers)
            
        # raise message if no primers were found for broad settings either
        primers = glob.glob(os.path.join(outprimers, ids +'.out'))
        if len(primers)==0:
            print("     No primers found for "+ids+" even with broader settings!")
            
        # raise message if there are errors in primer design
        error_file = os.path.join(outprimers, ids +'.err')
        file = open(error_file, 'r')
        errErrors = file.readlines()
        file.close()
        
        outErrors = []
        with open(os.path.join(outprimers, ids+'.out')) as file:
            for line in file.readlines():
                if "ERROR" in line:
                    outErrors.append(line)
                    
        if len(errErrors)>0 | len(outErrors)>0:
            print("     ERROR in primer design for "+ids+"!")
        
        # progress tracking
        if row%100 == 0:
            print("      primers designed for "+str(row)+" sequences")




if __name__=="__main__":
   # if len(sys.argv) != 3:
   #     print("Batch primer design takes three arguments: 1) the template CSV, 2) the output directory, and 3) filepath to primer3_core")
   #     print("1_primer3_batch_design.sh <TEMPLATES> <OUTDIR> <PRIMER3>")
   # else:
        main(sys.argv[1], sys.argv[2], sys.argv[3])
