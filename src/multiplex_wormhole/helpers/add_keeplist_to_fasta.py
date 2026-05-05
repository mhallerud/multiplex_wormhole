#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: ADD KEEPLIST TO FASTA
Purpose: Combines FASTA files for one record per unique primer, based on primer IDs.

Created on Thu Mar 14 15:00:12 2024
@author: maggiehallerud
"""

# load dependencies
import os
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
    # check inputs
    if not os.path.exists(MAIN_FA):
        raise InputError("MAIN_FA file could not be found!")
    if not os.path.exists(KEEPLIST_FA):
        raise InputError("KEEPLIST_FA file could not be found!")
    
    # set output filepath if not provided
    if OUTPATH==None:
    	OUTPATH= MAIN_FA.split('.fa')[0] + '_plusKeeplist.fa'

    # read in fastas
    main_primers, main_ids, main_seqids = readFasta(MAIN_FA)
    print(str(len(main_ids))+" MAIN sequences read.")
    keeplist_primers, keeplist_ids, keeplist_seqids = readFasta(KEEPLIST_FA)
    print(str(len(keeplist_ids))+" KEEPLIST sequences read.")
   
    # check for duplicate ids in each file
    if len(keeplist_ids) != len(set(keeplist_ids)):
        raise Exception("Non-unique primer IDs found in KEEPLIST_FA! Sequence headers in FASTA must be unique.")
    if len(main_ids) != len(set(main_ids)):
        raise Exception("Non-unique primer IDs found in MAIN_FA! Sequence headers in FASTA must be unique.")
    
    # if MAIN_FA contains pairIDs matching those in KEEPLIST_FA, remove from MAIN_FA to avoid duplicates
    # identify duplicate ids, sequences
    dup_ids = [main_seqids[x] for x in range(len(main_seqids)) if main_seqids[x] in keeplist_seqids]
    dup_seqs = [main_primers[x] for x in range(len(main_primers)) if main_primers[x] in keeplist_primers]
    if len(dup_seqs)>0:
        dup_seqs = [main_seqids[x] for x in range(len(main_seqids)) if main_primers[x] in dup_seqs]
        dup_ids.extend(dup_seqs)
    # remove primers belonging to templates in keeplist from the main FA
    if len(dup_ids)>0:
        dup_seqids = [x.split(".")[0] for x in dup_ids]
        main_primers = [main_primers[x] for x in range(len(main_primers)) if main_seqids[x] not in dup_seqids]
        main_ids = [main_ids[x] for x in range(len(main_ids)) if main_seqids[x] not in dup_seqids]
    
    # combine fastas and write to file
    combo_seqs = main_primers + keeplist_primers
    combo_ids = main_ids + keeplist_ids
    with open(OUTPATH, 'w') as file:
        for row in range(len(combo_seqs)):
            file.write(combo_ids[row])
            file.write(combo_seqs[row])
    print("     Output saved to: "+OUTPATH)



def readFasta(FA):
    # empty array to hold ids & seqs
    primerids = []    
    seqs = []
    # read in lines
    with open(FA, 'r') as file:
        lines = file.readlines()
        # make list of ids and seqs
        for line in lines:
            if '>' in line:
                primerids.append(line)
            else:
                seqs.append(line)
    # grab pair IDs
    seqid = [x.split(".")[0] for x in primerids]
    return seqs, primerids, seqid



class InputError(Exception):
    pass



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
