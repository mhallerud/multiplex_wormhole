#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: CONVERT CSV TO FASTA
Purpose: Converts input sequences from CSV to FASTA format

Created on Thu Dec 21 20:06:00 2023
@author: maggiehallerud
"""

# load modules
import sys
import os
import csv
import argparse
import traceback



def main(IN_CSV, OUT_FA, ID_FIELD="PrimerID", SEQ_FIELD="Sequence", 
         ENCODING=sys.getfilesystemencoding()):#TRY 'utf-8-sig' if fails
    """
    IN_CSV : CSV containing primer IDs and sequences [filepath]
    OUT_FA : FASTA file to be output [filepath]
    ID_FIELD : field name containing sequence IDs [string] default: PrimerID
    SEQ_FIELD : field name containing sequences [string] default: Sequence
    ENCODING : encoding of input CSV [string] default: python default
    -------
    Converts input sequences in CSV to FASTA format
    """
    if not os.path.exists(IN_CSV):
        raise Exception("IN_CSV file could not be found!")
    
    # check that output ends with FASTA or FA
    if not (OUT_FA.endswith(".fa") or OUT_FA.endswith(".fasta")):
        raise InputError("OUT_FA must end with .fa or .fasta extension")
   
    # set up empty arrays to hold sequences and ids separately
    sequences = []
    ids = []

    # read in csv
    try:
        with open(IN_CSV, 'r', encoding=ENCODING) as file:
            reader = csv.reader(file)
            # grab column index for each field
            header = next(reader)
            seqi = [header.index(i) for i in [SEQ_FIELD]][0]
            idi = [header.index(i) for i in [ID_FIELD]][0]
            for line in reader:
                # extract sequence and id based on field numbers
                sequences.append(line[seqi])
                ids.append(line[idi])
    except Exception:
        try:
            with open(IN_CSV, 'r', encoding='utf-8-sig') as file:
                reader = csv.reader(file)
                # grab column index for each field
                header = next(reader)
                seqi = [header.index(i) for i in [SEQ_FIELD]][0]
                idi = [header.index(i) for i in [ID_FIELD]][0]
                for line in reader:
                    # extract sequence and id based on field numbers
                    sequences.append(line[seqi])
                    ids.append(line[idi])
        except Exception as err:
            print(err)
            traceback.print_exc()
    
    # export data to fasta
    with open(OUT_FA, 'w') as file:
        for i in range(len(sequences)):
            file.write(">" + str(ids[i]) + "\n")
            file.write(sequences[i] + "\n")



# define InputError as an exception
class InputError(Exception):
    pass



def parse_args():
    # initialize argparser
    parser = argparse.ArgumentParser(description="Converts a CSV of primers to a FASTA. "+
                                     "Full documentation: https://mhallerud.github.io/multiplex_wormhole/csv-to-fasta",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # add required arguments
    parser.add_argument("-i", "--incsv", type=str, required=True)
    parser.add_argument("-o", "--outfa", type=str, required=True)
    # add optional arguments
    parser.add_argument("-p", "--primeridfield", type=str, default="PrimerID")
    parser.add_argument("-s", "--seqfield", type=str, default="Sequence")
    parser.add_argument("-e", "--csv_encoding", type=str, default=sys.getfilesystemencoding())
    return parser.parse_args()


         
def cli():
    # set up to run via command line
    args = parse_args()
    main(IN_CSV=args.incsv,
         OUT_FA=args.outfa, 
         ID_FIELD=args.primeridfield,
         SEQ_FIELD=args.seqfield,
         ENCODING=args.csv_encoding)



if __name__=="__main__":
    cli()
