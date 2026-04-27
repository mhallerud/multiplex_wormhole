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
    # set output filepath if not provided
    if OUTPATH==None:
    	OUTPATH= MAIN_FA.split('.fa')[0] + '_plusKeeplist.fa'

    # read in fastas
    main_seqs, main_ids = readFasta(MAIN_FA)
    print(str(len(main_ids))+" MAIN sequences read.")
    keeplist_seqs, keeplist_ids = readFasta(KEEPLIST_FA)
    print(str(len(keeplist_ids))+" KEEPLIST sequences read.")
   
    # check for duplicate ids in each file
    if len(keeplist_ids) != len(set(keeplist_ids)):
        raise Exception("Non-unique primer IDs found in KEEPLIST_FA! Sequence headers in FASTA must be unique.")
    if len(main_ids) != len(set(main_ids)):
        raise Exception("Non-unique primer IDs found in MAIN_FA! Sequence headers in FASTA must be unique.")
    
    # if MAIN_FA contains primerIDs matching those in KEEPLIST_FA, remove from MAIN_FA to avoid duplicates
    dup_ids = [keeplist_ids[x] for x in range(len(keeplist_ids)) if keeplist_ids[x] in main_ids]
    if len(dup_ids)>0:
        main_seqs = [main_seqs[x] for x in range(len(main_seqs)) if main_ids[x] not in dup_ids]
        main_ids = [main_ids[x] for x in range(len(main_ids)) if main_ids[x] not in dup_ids]
        
    
    # combine fastas and write to file
    combo_seqs = main_seqs + keeplist_seqs
    combo_ids = main_ids + keeplist_ids
    with open(OUTPATH, 'w') as file:
        for row in range(len(combo_seqs)):
            file.write(combo_ids[row])
            file.write(combo_seqs[row])
    print("     Output saved to: "+OUTPATH)


def readFasta(FA):
    # empty array to hold ids & seqs
    ids = []    
    seqs = []
    # read in lines
    with open(FA, 'r') as file:
        lines = file.readlines()
        # make list of ids and seqs
        for line in lines:
            if '>' in line:
                ids.append(line)
            else:
                seqs.append(line)
    return seqs, ids



if __name__=="__main__":
    main(sys.argv[1], sys.argv[2])
