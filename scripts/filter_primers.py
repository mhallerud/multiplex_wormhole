#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: FILTER PRIMERS
Created on Tue Dec 19 21:53:08 2023
@author: maggiehallerud

Purpose: Filters primer3 primer design output to remove primer pairs with likely secondary structures

Inputs: Path to primer directory, output path, and filtering parameters

Outputs: 1) Filtered primers CSV
         2) LocusIDs text file
"""

# load dependencies
import os
import sys
import glob
import csv
import re



def main(PRIMER_DIR, OUTPATH, Tm_LIMIT=45, dG_HAIRPINS=-2000, dG_END_LIMIT=-5000, dG_MID_LIMIT=-10000):
    """
    PRIMER_DIR : path to primer directory
        primer directory
    OUTPATH : string
        prefix for output files
    Tm_LIMIT : degrees Celsius
        minimum melting temperature allowed (default 45)
    dG_HAIRPINS : cal/mol
        maximum delta G allowed for hairpin structures (default -2000)
    dG_END_LIMIT : cal/mol
        maximium delta G allowed for dimers at ends of primers (default -5000)
    dG_MID_LIMIT : cal/mol
        maximum delta G allowed for dimers not at primer ends (default -10000)
    -------
    Removes primer pairs with predicted secondary structures;
    Outputs filtered primer pair details to a CSV
    """

    # make list of all primer3 output files
    locus_files = glob.glob(os.path.join(PRIMER_DIR,'*.out'))
    
    # raise error if no primer3 output files in PRIMER_DIR
    if len(locus_files)==0:
        raise InputError("PRIMER_DIR must contain .out and .err files from primer3 primer design")
            
    
    # set up an empty array to store primers that pass filtering
    filtered_primers=[]
    filtered_primers.append(["PrimerID","LocusID","PrimerPair","Direction","Sequence","StartBP","Length","AnnealingTempC","PropBound","AmpliconSize"])
    ids=[] # array for locus IDs that pass
    
    failed_loci = []
    
    for locus in locus_files:
        try:
            # get locus ID from pathname
            locusID = os.path.basename(locus).split('.')[0]
            
            # progress tracking for # primers per locus
            passed=0
            
            # read in file
            lines = []
            with open(locus, 'r') as file:
                for line in file.readlines():
                    lines.append(line.strip()) #.strip removes whitespace
            
            # find number of primer pairs
            primerpairs_line = list(filter(re.compile("PRIMER_PAIR_NUM_RETURNED*.").match, lines))[0]
            primerpairs = int(primerpairs_line.split("=")[1])
            
            # proceed if there are primers designed for this locus
            if primerpairs > 0:
                print("Testing "+locusID+"..............................")
                
                # loop through each primer pair to test criteria
                for N in range(primerpairs):
                    N = str(N) # convert to string
                    #print("     Primer pair "+N)
                    
                    # start counter for tests
                    tests=0
    
                    # check FW and REV primers separately for self-dimers
                    for DIR in ["LEFT","RIGHT"]:
                        # set strings based on direction
                        if DIR=="LEFT":
                            DIRNAME="FW"
                        elif DIR=="RIGHT":
                            DIRNAME="REV"
                        else:
                            pass
                        
                        # check for dimers with self
                        tests = checkStructure(lines, "PRIMER_"+DIR+"_"+N+"_SELF_ANY_STUCT", 
                                       Tm_LIMIT, dG_MID_LIMIT, DIRNAME+" Dimer with self (middle)", tests)                    
                        # check for dimers with self at ends
                        tests = checkStructure(lines, "PRIMER_"+DIR+"_"+N+"_SELF_END_STUCT", 
                                       Tm_LIMIT, dG_END_LIMIT, DIRNAME+" Dimer with self (end)", tests)
                        # check for hairpin dimers
                        tests = checkStructure(lines, "PRIMER_"+DIR+"_"+N+"_HAIRPIN_STUCT", 
                                       Tm_LIMIT, dG_HAIRPINS, DIRNAME+" Hairpin dimer", tests)
                    
                    # check for secondary structures between primers
                    tests = checkStructure(lines, "PRIMER_PAIR_"+N+"_COMPL_ANY_STUCT",
                                           Tm_LIMIT, dG_MID_LIMIT, DIRNAME+" Primer pair dimer (middle)", tests)                
                    # check for secondary structures between primer ends
                    tests = checkStructure(lines, "PRIMER_PAIR_"+N+"_COMPL_END_STUCT", 
                                           Tm_LIMIT, dG_END_LIMIT, DIRNAME+" Primer pair dimer (ends)", tests)
                    
                    # if all tests are passed (n=8), then the primer pairs are added to the output array
                    if tests==8:
                        # extract elements we want to save
                        FWseq_line = list(filter(re.compile("PRIMER_LEFT_"+N+"_SEQUENCE").match, lines))[0]
                        FWseq = FWseq_line.split("=")[1].lower()
                        FWpos_line = list(filter(re.compile("PRIMER_LEFT_"+N+"=").match, lines))[0]
                        FWpos = FWpos_line.split("=")[1].split(",")
                        FWtm_line = list(filter(re.compile("PRIMER_LEFT_"+N+"_TM").match, lines))[0]
                        FWtm = FWtm_line.split("=")[1]
                        FWbound_line = list(filter(re.compile("PRIMER_LEFT_"+N+"_BOUND").match, lines))[0]
                        FWbound = FWbound_line.split("=")[1]
                        REVseq_line = list(filter(re.compile("PRIMER_RIGHT_"+N+"_SEQUENCE").match, lines))[0]
                        REVseq = REVseq_line.split("=")[1].lower()
                        REVpos_line = list(filter(re.compile("PRIMER_RIGHT_"+N+"=").match, lines))[0]
                        REVpos = REVpos_line.split("=")[1]
                        REVtm_line = list(filter(re.compile("PRIMER_RIGHT_"+N+"_TM").match, lines))[0]
                        REVtm = REVtm_line.split("=")[1]
                        REVbound_line = list(filter(re.compile("PRIMER_RIGHT_"+N+"_BOUND").match, lines))[0]
                        REVbound = REVbound_line.split("=")[1]
                        ampSize_line = list(filter(re.compile("PRIMER_PAIR_"+N+"_PRODUCT_SIZE").match, lines))[0]
                        ampSize = ampSize_line.split("=")[1]
                    
                        # add FW and REV primers separately to the array
                        filtered_primers.append([locusID+"_"+N+"_FW", # primer ID
                                                 locusID, # locus ID
                                                 N, # primer pair #
                                                 "FW", # primer direction
                                                 FWseq, # primer sequence
                                                 FWpos[0], # primer position (start BP)
                                                 FWpos[1], # primer position (length)
                                                 FWtm, # primer annealing temp
                                                 FWbound, # primer proportion bound
                                                 ampSize]) # amplicon size
                        filtered_primers.append([locusID+"_"+N+"_REV", # primer ID
                                                 locusID, # locus ID
                                                 N, # primer pair #
                                                 "REV", # primer direction
                                                 REVseq, # primer sequence
                                                 REVpos[0], # primer position (start BP)
                                                 REVpos[1], # primer position (length)
                                                 REVtm, # primer annealing temp
                                                 REVbound, # primer proportion bound
                                                 ampSize]) # amplicon size
                        # progress tracking for the number of primer pairs per locus that pass
                        passed+=1
                        
            # add locus ID to list if any primer pairs passed
            if passed>0:
                ids.append(locusID)
            
            # clean up before restarting loop
            del lines
    
        except Exception:
            failed_loci.append(locusID)
    
    if len(failed_loci)>0:
        print("Filtering failed for the following loci. Retry manually if you want these included:")
        for l in failed_loci:
            print("          "+l)
    
    # Export filtered primers as CSV
    OUTCSV= OUTPATH + '.csv'
    if os.path.exists(OUTCSV):
        os.remove(OUTCSV) # remove if it already exists
    
    with open(OUTCSV, 'w', newline="\n") as file:
        writer = csv.writer(file)
        for row in filtered_primers:
            writer.writerow(row)
    
    # export locus IDs that passed as text file
    OUTIDS = OUTPATH+"_LocusIDs.txt"
    if os.path.exists(OUTIDS):
        os.remove(OUTIDS) # remove if it already exists
    
    with open(OUTIDS, 'w', newline="\n") as file:
        for row in ids:
            file.write(row+'\n')



def checkStructure(lines, structure_name, Tm_threshold, dG_threshold, message, tests):
    Structure = list(filter(re.compile(structure_name).match, lines))
    # if self secondary structures exist, test criteria
    if len(Structure) > 0:
        # extract annealing temp of structure
        Tm = Structure[0].split("=")[1].split(";")[0].split(":")[1].split("&")[0].strip()
        Tm_num = float(Tm)
        # extract deltaG of structure
        dG = Structure[0].split("=")[1].split(";")[1].split(" ")[3]
        dG_num = float(dG)
        # check if thresholds are met
        if Tm_num > Tm_threshold and dG_num < dG_threshold:
            #print("         "+message)
            #print("             "+structure_name+" failed tests")
            pass
        else:
            tests+=1
    else:
        tests+=1
    return tests



# define InputError as an exception
class InputError(Exception):
    pass



if __name__=="__main__":
    main(sys.argv[1],
         sys.argv[2],
         float(sys.argv[3]),
         float(sys.argv[4]),
         float(sys.argv[5]),
         float(sys.argv[6]))
