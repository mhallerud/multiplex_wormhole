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
import glob
import logging
import traceback
from datetime import datetime
import pandas
import numpy
import argparse
import subprocess


# load multiplex wormhole scripts
#sys.path.append('/Users/maggiehallerud/Desktop/multiplex_wormhole')#change to YOUR multiplex_wormhole path
sys.path.append(os.path.dirname(__file__))
from tabulate_dimers import main as tabulateDimers
from helpers.CSVtoFasta import main as csv2fasta
from helpers.logging_setup import setup_logging

## FIND PATH TO BINARY DEPENDENCIES
## NO SPACES ALLOWED IN PATHS- OTHERWISE CALLING FUNCTIONS WILL BREAK!
from helpers.setup_mfeprimer import main as setup_mfeprimer
MFEprimer_PATH = setup_mfeprimer()



def main(PRIMERS, ALL_DIMERS_dG=-8, END_DIMERS_dG=-4, BAD_DIMERS_dG=-10):
    """
    PRIMERS : FASTA or CSV (ID,Sequence) of primer sequences
    ALL_DIMERS_dG : DeltaG threshold used for dimer prediction (all dimers)
    END_DIMERS_dG : DeltaG threshold used for dimer prediction (end dimers)
    BAD_DIMERS_dG : DeltaG threshold used to count "bad" dimers
    -------------
    Calculates predicted dimer load and primer pairs involved, returns dimer output files.
    """    
    # setup logging
    logger = setup_logging(PRIMERS.split(".")[0]+".log", True, PRIMERS.split(".")[0])
    logger.info("START TIME: %s", datetime.now().strftime('%m/%d/%Y %I:%M:%S %p'))
    logger.info("")
    logger.info("panel_assessment inputs: ")
    logger.info("      PRIMERS: %s", PRIMERS)
    logger.info("      ALL_DIMERS_dG: %s", ALL_DIMERS_dG)
    logger.info("      END_DIMERS_dG: %s", END_DIMERS_dG)
    logger.info("      BAD_DIMERS_dG: %s", BAD_DIMERS_dG)
    logger.info("")

    # check paths
    if not os.path.exists(MFEprimer_PATH):
        raise InputError("Could not find MFEprimer!"+
                         "Set path on line 44 of panel_assessment.py which can be found here: "+
                         os.path.dirname(__file__))
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
    logger.info("Predicting dimers....")
    PREFIX = INPUT.replace(".fasta","").replace(".fa","")
    ALL_DIMERS = PREFIX+"_MFEdimers.txt"
    END_DIMERS = PREFIX+"_MFEdimers_ends.txt"
    try:
        subprocess.call(MFEprimer_PATH+" dimer -i "+INPUT+" -o "+ALL_DIMERS+" -d "+str(ALL_DIMERS_dG)+\
                        " -s 3 -m 50 --diva 3.8 --mono 50 --dntp 0.25 --oligo 50", shell=True)
        subprocess.call(MFEprimer_PATH+" dimer -i "+INPUT+" -o "+END_DIMERS+" -d "+str(END_DIMERS_dG)+\
                        " -s 3 -m 70 --diva 3.8 --mono 50 --dntp 0.25 --oligo 50 -p", shell=True)
    except Exception as err:
        logger.info("MFEprimer failed! Full error message:")
        logger.info(traceback.format_exc(err))
    
    # tabulate dimers into pairwise dimer table
    logger.info("Tabulating dimers (standard)...")
    OUT_DIMERS = PREFIX+"_PrimerPairDimers"
    tabulateDimers(ALL_DIMERS, 
                   END_DIMERS, 
                   OUT_DIMERS, 
                   "False",#specify if you want interactions per primer
                   deltaG=False)#count dimers
    logger.info("Tabulating dimers (mean min deltaG)")
    OUT_DELTAG = PREFIX+"_PrimerPairDeltaG"
    tabulateDimers(ALL_DIMERS,
                   END_DIMERS,
                   OUT_DELTAG,
                   "False",
                   deltaG=True)
    
    # count dimers
    logger.info("")
    logger.info("---------PANEL ASSESSMENT---------")
    countDimers(OUT_DIMERS+"_wide.csv", OUT_DELTAG+"_wide.csv", BAD_DIMERS_dG, logger)
    
    # close out logger
    logger.info("END TIME: %s", datetime.now().strftime('%m/%d/%Y %I:%M:%S %p'))
    logging.shutdown()



def countDimers(DIMERS_WIDE, DELTAG_WIDE, dG_THRESHOLD, logger): 
    ## SUMMARIZE DIMER COUNTS
    # read in pairwise dimer load CSV
    df = pandas.read_csv(DIMERS_WIDE)
    logger.info("Number of primer pairs assessed: %s", str(len(df)))
    # take sum across upper triangle, including diagonal
    s = numpy.triu(df,1).sum()
    logger.info("Number of pairwise dimers: %s", str(s))
    # take sum across rows (i.e., # interactions per primer pair)
    df = df.sum()[1:]#take sum across rows
    # sum of sums>0
    p = (df>0).sum()
    logger.info("Number of primer pairs involved in dimer interactions: %s", str(p))
    logger.info("Dimers per pair: "+str(round( s/len(df), 2)))
    
    ## SUMMARIZE DELTAG
    df = pandas.read_csv(DELTAG_WIDE)
    # take mean across upper triangle, including diagonal
    s = numpy.triu(df,1).mean()
    logger.info("Mean deltaG: "+str(round(s,2)))
    # convert to binary for primer pairs with dimers below a given deltaG threshold
    df.drop(df.columns[0], axis=1, inplace=True)
    df = df<dG_THRESHOLD
    df = df.replace({True: 1, False: 0})
    # take count of pairwise interactions with dimers below threshold
    b = numpy.triu(df,1).sum()
    logger.info("Number of pairwise interactions with 'bad' dimers: %s", str(b))
    logger.info("    (based on a deltaG threshold %s)", str(dG_THRESHOLD))



class InputError(Exception):
    pass



def parse_args():
    # initialize argparser
    parser = argparse.ArgumentParser()
    # add required arguments
    parser.add_argument("-i", "--input", type=str, required=True)
    # add optional arguments
    parser.add_argument("-a", "--alldimers_dg", type=float, default=-5)
    parser.add_argument("-e", "--enddimers_dg", type=float, default=-3)
    parser.add_argument("-b", "--baddimers_dg", type=float, default=-8)
    return parser.parse_args()



if __name__ == "__main__":
    # parse command-line arguments
    args = parse_args()
    # run panel assessment
    main(args.input,
         args.alldimers_dg,
         args.enddimers_dg,
         args.baddimers_dg)
