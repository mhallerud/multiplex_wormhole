#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: ADD KEEPLIST TO FASTA
Purpose: Combines FASTA files for one record per unique primer, based on primer IDs.

Created on Thu Mar 14 15:00:12 2024
@author: maggiehallerud
"""

# load dependencies
import argparse



def main(MAIN_FA, KEEPLIST_FA, OUTPATH=None):
    """
    Parameters
    ----------
    MAIN_FA : Template primer sequences [FASTA]
    KEEPLIST_FA : Keeplist primer sequences [FASTA]
    OUTPATH : Filepath to save merged FASTA [default: MAIN_FA+"_plusKeeplist.fa"]
    -------
    Outputs combined FASTA, attempts to handle redundant IDs / sequences
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
    
    # if MAIN_FA contains sequences matching those in KEEPLIST_FA, remove from MAIN_FA to avoid duplicates
    dup_seqs = [keeplist_seqs[x] for x in range(len(keeplist_seqs)) if keeplist_seqs[x] in main_seqs]
    if len(dup_seqs)>0:
        main_seqs = [main_seqs[x] for x in range(len(main_seqs)) if keeplist_seqs[x] in dup_seqs]
        main_ids = [main_ids[x] for x in range(len(main_ids)) if keeplist_seqs[x] in dup_seqs]
    
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



def parse_args():
    # initialize argparser
    parser = argparse.ArgumentParser()
    # add required arguments
    parser.add_argument("-i", "--input", type=str, required=True)
    parser.add_argument("-k", "--keeplist", type=str, required=True)
    # add optional arguments
    parser.add_argument("-o", "--outpath", type=str, default=None)
    return parser.parse_args()



if __name__=="__main__":
    args = parse_args()
    main(MAIN_FA=args.input, 
         KEEPLIST_FA=args.keeplist, 
         OUTPATH=args.outpath)
