#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: Add Keeplist Primers to Fasta
Input: Keeplist primers fasta and main primers fasta
Purpose: Combines FASTA files so that there is one record per unique primer, based on primer IDs.
Created on Thu Mar 14 15:00:12 2024

@author: maggiehallerud
"""

# load dependencies
import sys



def main(MAIN_FA, KEEPLIST_FA, OUTPATH=None):
    """
    Parameters
    ----------
    MAIN_FA : Fasta
        Specificity check passed output from step 3
    KEEPLIST_FA : Fasta
        Previously designed primers to add to newly designed primers

    Returns
    -------
    Outputs combined FASTA

    """
    # read in fastas
    main, main_ids = readFasta(MAIN_FA)
    keeplist, keeplist_ids = readFasta(KEEPLIST_FA)
   
    # check for duplicate ids
    all_ids = main_ids + keeplist_ids
    if len(keeplist_ids) != len(set(keeplist_ids)):
        raise Exception("Non-unique primer IDs found in KEEPLIST_FA! Fix IDs and then rerun.")
    elif len(set(all_ids)) != len(all_ids):
        raise Exception("KEEPLIST_FA has primer IDs that are also found in main fasta! Rename non-unique KEEPLIST primers IDs, then try again.")
    else:
        pass
        
    # combine fastas and write to file
    combo_fa = main + keeplist
    if OUTPATH=None:
    	OUTPATH= MAIN_FA.split('.')[0] + '_plusKeeplist.fa'
    with open(OUTPATH, 'w') as file:
        for line in combo_fa:
            file.write(line)


def readFasta(FA):
    # empty array to hold ids
    ids = []    
    # read in lines
    with open(FA, 'r') as file:
        lines = file.readlines()
        # make list of ids
        for line in lines:
            if '>' in line:
                ids.append(line)
    return lines, ids



if __name__=="__main__":
    main(sys.argv[1], sys.argv[2])
