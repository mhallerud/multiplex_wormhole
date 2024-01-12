#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: CHECK PRIMER SPECIFICITY
Created on Wed Dec 20 00:14:22 2023
@author: maggiehallerud

Purpose: Check primer specificity against all other available loci and/or genome files
        This ensures that passing primers will only amplify one target locus
Inputs: 1) Primer CSV output from filtered_primers.py
        2) Template file to check specificity against
        3) Prefix for output filepaths
Outputs: 1) CSV with all primers that passed the specificity check
         2) FASTA with all primers that passed the specificity check
         3) CSV file with primer info for all primers that failed the check
"""

# load dependencies
import os
import sys
import csv



def main(PRIMERS, TARGET, OUTPATH):
    """
    PRIMERS : CSV
        Primer details from filter_primers.py
    TARGET : CSV or FASTQ/FASTA/FASTQ.GZ
        Sequences from DNA available to primers for amplification
    OUTPATH : filepath
        path and prefix for outputs
    -------
    Checks specificity of PRIMERS against TARGET
    Returns CSVs with primers passing and failing the specificity check
    """
      
    # read in primer CSV
    with open(PRIMERS, 'r', newline='\n') as file:
        reader = csv.reader(file)
        next(reader) # skip header
        primers = list(reader)
    
    # empty arrays to hold passed and failed primer info
    passed = []
    passed_seq = []
    failed = []
    passed.append(["PrimerID","LocusID","PrimerPair","Direction","Sequence","StartBP","Length","AnnealingTempC","PropBound","AmpliconSize"])
    failed.append(["PrimerID","LocusID","PrimerPair","Direction","Sequence","StartBP","Length","AnnealingTempC","PropBound","AmpliconSize"])
    
    # use zgrep if input is a *.fastq.gz file
    if TARGET.endswith("fastq.gz") or TARGET.endswith("fq.gz"):
        # check for each primer
        # check for each primer
        for i in range(len(primers)):
            primer_seq = primers[i][4]# get primer sequence
            # remove adapter sequence and make uppercase
            primer_seq = primer_seq.replace('tcgtcggcagcgtcagatgtgtataagagacag','')
            primer_seq = primer_seq.replace('gtctcgtgggctcggagatgtgtataagagacag','')
            primer_seq = primer_seq.upper()
            # search for sequence in file using grep, save to temp file
            if os.path.exists("grep.tmp"):
                os.remove("grep.tmp")
            zgrep = "zgrep "+primer_seq+" "+TARGET+" > grep.tmp"
            os.system(zgrep)
            # read in matching sequences
            matches = [i for i in open("grep.tmp").readlines()]
            # if there is more than one match, the primer fails
            if len(matches)>1:
                failed.append(primers[i])
            # if there is 1 or 0 matches, primer passes
            # NOTE: consider changing this to require 1 match
            if len(matches)<=1:
                passed.append(primers[i])
                primerID = primers[i][0]
                passed_seq.append(">"+primerID)
                passed_seq.append(primer_seq)
    # use normal grep if input is a *.fasta or *.fastq file
    elif TARGET.endswith(".fasta") or TARGET.endswith(".fastq") or TARGET.endswith(".fa") or TARGET.endswith(".fq") or TARGET.endswith(".csv"):
        # check for each primer
        for i in range(len(primers)):
            primer_seq = primers[i][4]# get primer sequence
            # remove adapter sequence and make uppercase
            primer_seq = primer_seq.replace('tcgtcggcagcgtcagatgtgtataagagacag','')
            primer_seq = primer_seq.replace('gtctcgtgggctcggagatgtgtataagagacag','')
            primer_seq = primer_seq.upper()
            # search for sequence in file using grep, save to temp file
            if os.path.exists("grep.tmp"):
                os.remove("grep.tmp")
            zgrep = "grep "+primer_seq+" "+TARGET+" > grep.tmp"
            os.system(zgrep)
            # read in matching sequences
            matches = [i for i in open("grep.tmp").readlines()]
            # if there is more than one match, the primer fails
            if len(matches)>1:
                failed.append(primers[i])
            # if there is 1 or 0 matches, primer passes
            # NOTE: consider changing this to require 1 match
            if len(matches)<=1:
                passed.append(primers[i])
                primerID = primers[i][0]
                passed_seq.append(">"+primerID)
                passed_seq.append(primer_seq)
    # if neither fastq.gz, fasta, or fastq file, raise error
    else:
        raise InputError("TARGET file of must be a CSV file (.csv) or a FASTA (.fasta, .fa), FASTQ (.fastq, .fq), or FASTQ.GZ file (.fastq.gz or .fq.gz)")
        
    # write passed primers to a CSV file
    passed_csv_path = OUTPATH+"_passed.csv"
    with open(passed_csv_path, 'w') as file:
        writer= csv.writer(file)
        for row in passed:
            writer.writerow(row)
    
    # write passed primer sequences to a FASTA file
    passed_fa_path = OUTPATH+"_passed.fa"
    with open(passed_fa_path, 'w') as file:
        for row in passed_seq:
            file.write(row+"\n")
    
    # write failed primers to CSV file
    failed_csv_path = OUTPATH+"_failed.csv"
    with open(failed_csv_path, 'w') as file:
        writer = csv.writer(file)
        for row in failed:
            writer.writerow(row)



def reverse_complement(sequence):
    tab = str.maketrans("acgt", "tgca")
    revcomp = sequence.translate(tab)[::-1]
    return revcomp



# define InputError as an exception
class InputError(Exception):
    pass



if __name__=="__main__":
    main(sys.argv[1],
         sys.argv[2],
         sys.argv[3])