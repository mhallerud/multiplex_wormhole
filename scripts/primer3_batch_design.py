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
import pandas as pd
import warnings



def main(IN_CSV, OUTDIR, PRIMER3_PATH):
    """
    IN_CSV : CSV containing template sequences in primer3 format with header
        field1=ID, field2=DNA sequence, field3=amplicon (start_bp,length)

    OUTDIR : directory path
        an output directory where primer details are saved
    
    PRIMER3_PATH : primer3_core path
    	location of primer3_core software
    ------
    Designs primers for all loci in IN_CSV; 
    Returns primer output files in OUTDIR/1_InitialPrimers
    """
    # check inputs before proceeding
    if not os.path.exists(IN_CSV):
        raise InputPathError("Could not find IN_CSV <"+IN_CSV+">")
    if not os.path.exists(PRIMER3_PATH):
        raise InputPathError("Could not find PRIMER3_PATH <"+PRIMER3_PATH+">")
    if not os.path.exists(OUTDIR):
        raise InputPathError("Could not find OUTDIR <"+OUTDIR+">")
    
    # find script directory 
    SCRIPTPATH=os.path.dirname(__file__)

    # define paths to primer3 settings
    primer3_sh = os.path.join(SCRIPTPATH, 'primer3.sh')
    basedir = os.path.dirname(SCRIPTPATH)
    strict = os.path.join(basedir, 'primer3_settings/primer3_Base_NoSecondaryFilters.txt')
    broad = os.path.join(basedir, 'primer3_settings/primer3_Broad_NoSecondaryFilters.txt')
    if not os.path.exists(strict):
        raise InputPathError("Could not find primer3 settings <"+strict+">." + \
                             "Did you move scripts or primer3_settings out of your multiplex_wormhole directory?")
    if not os.path.exists(broad):
        raise InputPathError("Could not find primer3 settings <"+broad+">." + \
                             "Did you move scripts or primer3_settings out of your multiplex_wormhole directory?")

    # read in input sequences
    print("Reading in sequences......")
    # check that input sequences follow the proper formatting
    ids, seqs, snps = readCheckInputCSV(IN_CSV)

    # set up output folder structure
    if not os.path.exists(OUTDIR):
        os.mkdir(OUTDIR)
    
    # set up output directory to save primers
    outprimers = os.path.join(OUTDIR,'1_InitialPrimers')
    # if the folder already exists, delete and remake
    # this avoids issues building up between runs
    if os.path.exists(outprimers):
        print("Removing previous outputs in "+outprimers+"...")
        shutil.rmtree(outprimers)
    os.mkdir(outprimers)
    
    # design FW and REV primers for each sequence
    print("Designing initial primers......")
    for row in range(len(ids)):        
        id = ids[row]
        # run primer3 based on strict settings
        #print(" - Designing primers for "  + ids)
        os.system(primer3_sh +' '+ PRIMER3_PATH +' '+ id +' '+ seqs[row] +' '+ snps[row] +' '+ strict +' '+ outprimers)
    
    	# rerun with broad settings if no primers were found
        primers = glob.glob(os.path.join(outprimers, id +'.out'))
        if len(primers)==0:
            print("     No primers found for "+ids[row]+"! Retrying with broader settings.")
            os.system(primer3_sh +' '+ PRIMER3_PATH + ' '+ id +' '+ seqs[row] +' '+ snps[row] +' '+ broad +' '+ outprimers)
            
        # raise message if no primers were found for broad settings either
        primers = glob.glob(os.path.join(outprimers, id +'.out'))
        if len(primers)==0:
            print("     No primers found for "+id+" even with broader settings!")
            
        # raise message if there are errors in primer design
        error_file = os.path.join(outprimers, id +'.err')
        file = open(error_file, 'r')
        errErrors = file.readlines()
        file.close()
        
        outErrors = []
        with open(os.path.join(outprimers, id+'.out')) as file:
            for line in file.readlines():
                if "ERROR" in line:
                    outErrors.append(line)
                    
        if len(errErrors)>0 | len(outErrors)>0:
            print("     ERROR in primer design for "+ids[row]+"!")
        
        # progress tracking
        if row%100 == 0:
            print("      primers designed for "+str(row)+" sequences")



# define class expected for each row in IN_CSV
def readCheckInputCSV(IN_CSV):
    try:
        # Read the CSV file using the specified delimiter and header settings
        df = pd.read_csv(IN_CSV, delimiter=",", header=0)

        # Check header names
        head = df.columns
        if 'SEQUENCE_ID' not in head:
            raise InputParseError("SEQUENCE_ID field missing from "+IN_CSV)
        if 'SEQUENCE_TEMPLATE' not in head:
            raise InputParseError("SEQUENCE_TEMPLATE field missing from "+IN_CSV)
        if 'SEQUENCE_TARGET' not in head:
            raise InputParseError("SEQUENCE_TARGET field missing from "+IN_CSV)
        
        # Check that sequence IDs don't have "."
        ids = df["SEQUENCE_ID"]
        if any(['.' in ids[i] for i in range(len(ids))]):
            raise InputParseError("Some SEQUENCE_ID values contain a period '.' - remove any periods from SEQUENCE_ID and try again.")
        # Check that sequence IDs are unique
        if len(set(ids))<len(ids):
            raise InputParseError("SEQUENCE_ID field contains non-unique values - rename or remove duplicates and try again.")

        # Check that templates only contain ACTGN characters
        seqs = df["SEQUENCE_TEMPLATE"]
        allowed = {'a','c','g','t','n','A','C','G','T','N'}
        if any([len(set(seqs[i])-allowed)!=0 for i in range(len(seqs))]):
            raise InputParseError("SEQUENCE_TEMPLATE field contains characters other than ACTGN or actgn - remove these other characters and try again.")
        
        # Check that targets follow the <start>,<length> format
        snps = df["SEQUENCE_TARGET"]
        if any([len(snps[i].split(","))!=2 for i in range(len(snps))]):
            raise InputParseError("Some SEQUENCE_TARGET values do not match the format STARTBP,LENGTH - fix and try again")
        
        # Check that targets don't extend past end of template sequence
        seqlens = [len(seqs[i]) for i in range(len(seqs))]
        startbps = [int(snps[i].split(",")[0]) for i in range(len(snps))]
        if any([startbps[i]>seqlens[i] for i in range(len(seqs))]):
            warnings.warn("Some target's startBP are beyond the sequence - no primers will be created." + \
                          "Check that your SEQUENCE_TARGET format is in the format STARTBP,LENGTH "+ \
                           "and that targets are paired with their corresponding sequences.", InputWarning)
        targlens = [int(snps[i].split(",")[1]) for i in range(len(snps))]
        endbps = [startbps[i]+targlens[i] for i in range(len(snps))]
        if any([endbps[i]>seqlens[i] for i in range(len(seqs))]):
            warnings.warn("Some target's length extends beyond sequence - no primers will be created." + \
                          "Check that your SEQUENCE_TARGET format is in the format STARTBP,LENGTH "+ \
                          "and that targets are paired with their corresponding sequences.", InputWarning)
        
        return([ids, seqs, snps])  # Return fields
    except pd.errors.ParserError:
        raise InputParseError("IN_CSV <"+IN_CSV+"> was found but could not be read!")




# define custom exception and warning types
class InputPathError(Exception):
    pass
class InputParseError(Exception):
    pass
class InputWarning(UserWarning):
    pass




if __name__=="__main__":
   # if len(sys.argv) != 3:
   #     print("Batch primer design takes three arguments: 1) the template CSV, 2) the output directory, and 3) filepath to primer3_core")
   #     print("1_primer3_batch_design.sh <TEMPLATES> <OUTDIR> <PRIMER3>")
   # else:
        main(sys.argv[1], sys.argv[2], sys.argv[3])
