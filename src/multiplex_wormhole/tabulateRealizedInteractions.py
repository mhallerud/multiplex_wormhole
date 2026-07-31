#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CONVERT GT-SEQ DIMERS TO TABLE

Created on Thu Jul 30 20:51:16 2026

Code for parsing and pivoting the primer interaction data was developed with the 
assistance of Claude (Anthropic, Sonnet 4.5, 2026).

@author: maggiehallerud
"""

# load dependencies
import os
import glob
import pandas as pd
import argparse



GTSEQ_DIMERS = "lane1-s865-index--CGCTCAGTTC-TCGTGGAGCG-PEPE97_Hoopa298_a_S865.assembled.fastq.hash.interactions.tsv"

def main(GTSEQ_DIMERS, OUTPREFIX, FILELIST=False):
    if not FILELIST:
        # read in file
        if not os.path.exists(GTSEQ_DIMERS):
            raise InputError(GTSEQ_DIMERS+" file could not be found!")
        with open(GTSEQ_DIMERS, 'r') as file:
            lines = file.readlines()
    else:
        # read in multiple files
        files = glob.glob(GTSEQ_DIMERS)
        if len(files)<1:
            raise InputError("Files not found at: "+GTSEQ_DIMERS)
        lines = []
        for f in files:
            with open(f, 'r'):
                lines.append(file.readlines())
    
    # find index for each component (multiples possible- esp. with mult. files)
    idx1 = [x for x in range(len(lines)) if lines[x]=="***Proper On-Target Primer combinations\n"]
    idx2 = [x for x in range(len(lines)) if lines[x]=="***Double-Complement Primer-Artifacts***\n"]
    idx3 = [x for x in range(len(lines)) if lines[x]=="****FWD-REV Primer mis-primes****\n"]
    idx4 = [x for x in range(len(lines)) if lines[x]=="******FWD-FWD Primer mis-primes********\n"]
    idx5 = [x for x in range(len(lines)) if lines[x]=="****REV-FWD Primer mis-primes****\n"]
    idx6 = [x for x in range(len(lines)) if lines[x]=="******REV-REV Primer mis-primes********\n"]
    
    # set up arrays to hold empties
    ontarget = []
    doublecompl = []
    fwdrevmisprime = []
    fwdfwdmisprime = []
    revfwdmisprime = []
    revrevmisprime = []
    
    # mispriming defns based on GT-seq:
    # on-target: matches probe
    # double complement: primer dimer- primers anneal directly to one another
    # mis-primes: i.e., off-target amplification
    for i in range(len(idx1)):
        ontarget.extend(lines[(idx1[i]+2) : idx2[i]])
        doublecompl.extend(lines[(idx2[i]+1) : idx3[i]])
        fwdrevmisprime.extend(lines[(idx3[i]+1) : idx4[i]])
        fwdfwdmisprime.extend(lines[(idx4[i]+1) : idx5[i]])
        revfwdmisprime.extend(lines[(idx5[i]+1) : idx6[i]])
        if i<(len(idx1)-1):
            revrevmisprime.extend(lines[(idx6[i]+1) : idx1[i+1]])
        else:
            revrevmisprime.extend(lines[(idx6[i]+1) : len(lines)])
    
    # combine misprimes
    misprimes = fwdrevmisprime + fwdfwdmisprime + revfwdmisprime + revrevmisprime
    
    # build the three matrices
    ontarget_matrix = pivotToMatrix(ontarget)
    doublecompl_matrix = pivotToMatrix(doublecompl)
    misprimes_matrix = pivotToMatrix(misprimes)
    
    # write each to its own CSV
    ontarget_matrix.to_csv(OUTPREFIX + "_ontarget_matrix.csv")
    doublecompl_matrix.to_csv(OUTPREFIX + "_doublecomplement_matrix.csv")
    misprimes_matrix.to_csv(OUTPREFIX + "_misprimes_matrix.csv")
    
    # log outputs
    print("On-target read count:", ontarget_matrix.sum().sum())
    print("Double-complement primer read count:", doublecompl_matrix.sum().sum().sum,)
    print("Mis-priming read count:", misprimes_matrix.sum().sum())



def pivotToMatrix(section_lines):
    """
    Convert a list of raw lines into an n x n primer1 (rows) x primer2 (cols)
    matrix of counts. Missing combinations are filled with 0.
    Duplicate primer1/primer2 pairs (if any) are summed.
    """
    df = parseLines(section_lines)
    if df.empty:
        return pd.DataFrame()
    matrix = df.pivot_table(
        index="primer1",
        columns="primer2",
        values="count",
        aggfunc="sum",
        fill_value=0,
    )
    return matrix



def parseLines(section_lines):
    """
    Parse lines of the form 'primer1 primer2\tcount' into a long-format
    DataFrame with columns primer1, primer2, count.
    """
    records = []
    for line in section_lines:
        line = line.strip()
        if not line:
            continue
        combo, count = line.split("\t")
        primer1, primer2 = combo.split(" ", 1)  # split on first space only
        records.append((primer1, primer2, int(count)))
    return pd.DataFrame(records, columns=["primer1", "primer2", "count"])



class InputError(Exception):
    pass


def parse_args():
    # initialize argparser
    parser = argparse.ArgumentParser(description="Converts GT-seq interactions output into a table. ",
                                     #"Full documentation: https://mhallerud.github.io/multiplex_wormhole/tabulate-dimers",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # add required arguments
    parser.add_argument("-i", "--gtseq-dimers", type=str, required=True,
                        help="Filepath to GTseq-Primer-Interaction-Test output files (use * for multiple files)")
    parser.add_argument("-o", "--outpath", type=str, required=True,
                        help="Prefix for output tables (including directory)")
    parser.add_argument("-l", "--filelist", action="store_true",
                        help="Input is for list of files?")
    return parser.parse_args()



def cli():
    # parse command line arguments
    args = parse_args()
    # run main
    main(GTSEQ_DIMERS = args.gtseq_dimers,
         OUTPREFIX = args.outpath,
         FILELIST = args.filelist)



if __name__=="__main__":
    cli()