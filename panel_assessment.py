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
from scripts.extras.CSVtoFasta import main as csv2fasta




def main(PRIMERS, ALL_DIMERS_G=-5, END_DIMERS_G=-3, BAD_DIMERS_G=-5):
    """
    PRIMERS : FASTA or CSV (ID,Sequence) of primer sequences
    ALL_DIMERS_G : DeltaG threshold used for dimer prediction (all dimers)
    END_DIMERS_G : DeltaG threshold used for dimer prediction (end dimers)
    BAD_DIMERS_G : DeltaG threshold used to count "bad" dimers
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
    print("Predicting dimers.....")
    PREFIX = INPUT.replace(".fasta","").replace(".fa","")
    ALL_DIMERS = PREFIX+"_MFEdimers.txt"
    END_DIMERS = PREFIX+"_MFEdimers_ends.txt"
    # lower these deltaG thresholds since we want to include "minor" dimers here too
    os.system(MFEprimer_PATH+" dimer -i "+INPUT+" -o "+ALL_DIMERS+" -d "+str(ALL_DIMERS_G)+
              " -s 3 -m 50 --diva 3.8 "+
              "--mono 50 --dntp 0.25 --oligo 50")
    os.system(MFEprimer_PATH+" dimer -i "+INPUT+" -o "+END_DIMERS+" -d "+str(END_DIMERS_G)+
              " -s 3 -m 70 --diva 3.8 "+
              "--mono 50 --dntp 0.25 --oligo 50 -p")
    
    # tabulate dimers into pairwise dimer table
    print("Tabulating dimers.......")
    tabulateDimers(ALL_DIMERS, 
                   END_DIMERS, 
                   PREFIX+"PrimerPairInteractions", 
                   "False")#specify if you want interactions per primer
    
    # count dimers
    print("")
    print("----------PANEL ASSESSMENT---------")
    countDimers(PREFIX+"PrimerPairInteractions_wide.csv", BAD_DIMERS_G)



def countDimers(DIMERS_WIDE, DELTA_G): 
    """
    Parameters
    ----------
    DIMERS_WIDE : CSV filepath of max deltaG for pairwise primer interactions
    DELTA_G : Upper DeltaG threshold for counting dimers (default: -5)
    -------
    Returns mean deltaG and predicted dimer load across panel.
    """
    # read in pairwise dimer load CSV
    df = pandas.read_csv(DIMERS_WIDE)
    print("Number of primer pairs assessed: "+str(len(df)))
    # take mean across upper triangle, including diagonal
    s = numpy.triu(df,1).mean()
    print("Mean deltaG: "+str(s))
    # convert to binary for primer pairs with dimers below a given deltaG threshold
    df.drop(df.columns[0], axis=1, inplace=True)
    df = df<DELTA_G
    df = df.replace({True: 1, False: 0})
    # take count of pairwise interactions with dimers below threshold
    b = numpy.triu(df,1).sum()
    print("Number of pairwise interactions with 'bad' dimers: "+str(b))
    print("    (based on a deltaG threshold "+str(DELTA_G)+")")
    # take sum across rows (i.e., # interactions per primer pair)
    df = df.sum()[1:]#take sum across rows
    # sum of sums<0
    p = (df>0).sum()
    print("Number of primer pairs involved in 'bad' dimer interactions: "+str(p))



class InputError(Exception):
    pass



if __name__=="__main__":
    main(sys.argv[0])