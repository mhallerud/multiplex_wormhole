#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: Add Whitelist Primers to Fasta
Input: whitelist primers fasta and main primers fasta
Created on Thu Mar 14 15:00:12 2024

@author: maggiehallerud
"""

# load dependencies
import sys



def main(MAIN_FA, WHITELIST_FA):
    """
    Parameters
    ----------
    MAIN_FA : Fasta
        Specificity check passed output from step 3
    WHITELIST_FA : Fasta
        Existing primer set to include

    Returns
    -------
    Output fasta in same folder as MAIN_FA

    """
    # read in fastas
    main, main_ids = readFasta(MAIN_FA)
    whitelist, whitelist_ids = readFasta(WHITELIST_FA)
   
    # check for duplicate ids
    all_ids = main_ids + whitelist_ids
    if len(whitelist_ids) != len(set(whitelist_ids)):
        raise Exception("Non-unique primer IDs found in whitelist fasta! Fix IDs and then rerun.")
    elif len(set(all_ids)) != len(all_ids):
        raise Exception("Whitelist has primer IDs that are also found in main fasta! Rename non-unique whitelist primers IDs, then try again.")
    else:
        pass
        
    # combine fastas and write to file
    combo_fa = main + whitelist
    outpath = MAIN_FA.split('.')[0] + '_plusWhitelist.fa'
    with open(outpath, 'w') as file:
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
