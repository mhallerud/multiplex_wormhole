#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: PANEL ASSESSMENT with multiplex wormhole

Created on Wed Apr  8 21:57:37 2026
@author: maggiehallerud

IMPORTANT NOTE: Before running this script, make sure that you have installed the 
dependency MFEprimer (available at https://www.mfeprimer.com/mfeprimer-3.1/#2-command-line-version)


Input preparation:
    PRIMERS: Input FASTA or CSV file containing primer sequences. 
        Sequence identifiers must include a unique primer pair ID + 
        the primer direction (.FWD/.REV or .FW/.RV) separated by a period,
        with no other periods in the header. 
        E.g. "MACA1.FW" and "MACA1.RV"
        
        If CSV format, must contain 'PrimerID' and 'Sequence' fields.
"""

# load packages
import os
import sys
import pandas
import numpy


## SET PATHS TO DEPENDENCIES:
MFEprimer_PATH='/Users/maggiehallerud/Marten_Primer_Design/Plate1_First55Pairs_Sep2023/mfeprimer-3.2.7-darwin-10.6-amd64'#full path to mfeprimer location

# load multiplex wormhole scripts
sys.path.append('/Users/maggiehallerud/Desktop/multiplex_wormhole')#change to YOUR multiplex_wormhole path
from scripts.tabulate_MFEprimer_dimers import main as tabulateDimers
from scripts.CSVtoFasta import main as csv2fasta




def main(PRIMERS):
    """
    PRIMERS : FASTA or CSV (ID,Sequence) of primer sequences
    -------------
    Calculates predicted dimer load and primer pairs involved, returns dimer output files.
    """
    # convert CSV primers to FASTA, if needed
    if(PRIMERS.endswith(".csv")):
        # check for requisite fields first...
        CSV = pandas.read_csv(PRIMERS)
        if 'PrimerID' not in list(CSV.columns) or 'Sequence' not in list(CSV.columns):
            raise InputError("CSV is missing PrimerID or Sequence field- fix and try again.")
        # convert CSV to FASTA
        FASTA = PRIMERS.replace(".csv", ".fasta")
        csv2fasta(PRIMERS, OUT_FA=FASTA, ID_FIELD="PrimerID", SEQ_FIELD="Sequence")
        # set input filepath for next step
        INPUT = FASTA
    else:
        INPUT = PRIMERS

    # predict dimers with MFEprimer dimer function
    PREFIX = INPUT.replace(".fasta","").replace(".fa","")
    ALL_DIMERS = PREFIX+"_MFEdimers.txt"
    END_DIMERS = PREFIX+"_MFEdimers_ends.txt"
    os.system(MFEprimer_PATH+" dimer -i "+INPUT+" -o "+ALL_DIMERS+" -d -8 -s 3 -m 50 --diva 3.8 "+
              "--mono 50 --dntp 0.25 --oligo 50")
    os.system(MFEprimer_PATH+" dimer -i "+INPUT+" -o "+END_DIMERS+" -d -5 -s 3 -m 70 --diva 3.8 "+
              "--mono 50 --dntp 0.25 --oligo 50 -p")
    
    # tabulate dimers into pairwise dimer table
    tabulateDimers(ALL_DIMERS, 
                   END_DIMERS, 
                   PREFIX+"PrimerPairInteractions", 
                   "False")#specify if you want interactions per primer
    
    # count dimers
    print("")
    countDimers(PREFIX+"PrimerPairInteractions_wide.csv")



def countDimers(DIMERS_WIDE): 
    # read in pairwise dimer load CSV
    df = pandas.read_csv(DIMERS_WIDE)
    print("Number of primer pairs assessed: "+str(len(df)))
    # take sum across upper triangle, including diagonal
    s = numpy.triu(df,1).sum()
    print("Number of pairwise dimers: "+str(s))
    # take sum across rows (i.e., # interactions per primer pair)
    df = df.sum()[1:]#take sum across rows
    # sum of sums>0
    p = (df>0).sum()
    print("Number of primer pairs involved in dimer interactions: "+str(p))



class InputError(Exception):
    pass



if __name__=="__main__":
    main(sys.argv[0])