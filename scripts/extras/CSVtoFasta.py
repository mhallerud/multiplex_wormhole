#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: Convert CSV to FASTA
Created on Thu Dec 21 20:06:00 2023
@author: maggiehallerud
"""

# load modules
import sys
import csv



def main(IN_CSV, OUT_FA, ID_FIELD="PrimerID", SEQ_FIELD="Sequence"):
    """
    IN_CSV : CSV containing primer IDs and sequences [filepath]
    OUT_FA : FASTA file to be output [filepath]
    ID_FIELD : field name containing sequence IDs [string] default: PrimerID
    SEQ_FIELD : field name containing sequences [string] default: Sequence
    -------
    Converts input sequences in CSV to FASTA format
    """
    
    
    # check that output ends with FASTA or FA
    if not (OUT_FA.endswith(".fa") or OUT_FA.endswith(".fasta")):
        raise InputError("OUT_FA must end with .fa or .fasta extension")
   
    # set up empty arrays to hold sequences and ids separately
    sequences = []
    ids = []

    # read in csv
    with open(IN_CSV, 'r') as file:
        reader = csv.reader(file)
        # grab column index for each field
        header = next(reader)
        seqi = [header.index(i) for i in [SEQ_FIELD]][0]
        idi = [header.index(i) for i in [ID_FIELD]][0]
        for line in reader:
            # extract sequence and id based on field numbers
            sequences.append(line[seqi])
            ids.append(line[idi])
    
    # export data to fasta
    with open(OUT_FA, 'w') as file:
        for i in range(len(sequences)):
            file.write(">" + str(ids[i]) + "\n")
            file.write(sequences[i] + "\n")



# define InputError as an exception
class InputError(Exception):
    pass



# set up to run via command line
if __name__=="__main__":
    main(sys.argv[1],
         sys.argv[2],
         sys.argv[3],
         sys.argv[4])