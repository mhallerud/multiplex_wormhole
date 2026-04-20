#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 14:31:28 2024

@author: maggiehallerud
"""

import traceback
import signal
import os
import sys
import glob
import shutil

sys.path.append(os.path.dirname(__file__))
from optimize_primers import main as optimizeMultiplex



def multipleOptimizations(N_RUNS, PRIMER_FA, DIMER_SUMS, DIMER_TABLE, OUTPATH, N_LOCI, 
                          KEEPLIST=None, TIMEOUT=5, VERBOSE=False, SEED=None,
                          SIMPLE=5000, ITERATIONS=10000, BURNIN=100, DECAY_RATE=0.98, 
                          T_INIT=None, T_FINAL=None, PARTITIONS=1000, DIMER_ADJ=0.1, PROB_ADJ=2):
    # set up empty array to hold overall dimer load 
    loads = [['Run', 'TotalDimers']]
    
    # while instead of for loop allows runs that timeout to restart
    run = 1
    while run <= N_RUNS:
        signal.alarm(TIMEOUT*60) # set timer- this helps to timeout runs with infinite loops in the optimization process
        try:
            try:
                cost = optimizeMultiplex(PRIMER_FASTA=PRIMER_FA,
                                         DIMER_SUMS=DIMER_SUMS, 
                                         DIMER_TABLE=DIMER_TABLE,
                                         OUTPATH=OUTPATH+"_Run"+str(run).zfill(2),
                                         N_LOCI=N_LOCI, 
                                         KEEPLIST=KEEPLIST,
                                         SEED=None,
                                         VERBOSE=False,
                                         SIMPLE=SIMPLE, 
                                         ITERATIONS=ITERATIONS,
                                         BURNIN=BURNIN,
                                         DECAY_RATE=DECAY_RATE,
                                         T_INIT=T_INIT,
                                         T_FINAL=T_FINAL,
                                         PARTITIONS=PARTITIONS,
                                         DIMER_ADJ=DIMER_ADJ,
                                         PROB_ADJ=PROB_ADJ,
                                         MAKEPLOT=False)
                loads.append([str(run), str(cost)])
            except Exception as e:
                print("AN ERROR OCCURRED IN RUN "+str(run))
                print(e)
                traceback.print_exc()
                print("---")
                run+=1
                continue
        except TimeoutException:
            print("TimeoutException!")
            continue
        # reset alarm
        else:
            signal.alarm(0)
        print(" ")
        run+=1 #move to next run (only if no timeout)
    
    # export dimer loads per run
    with open(OUTPATH+'_RunSummary.csv', 'w') as file:
        for row in loads:
            file.write(row[0]+','+row[1]+"\n")
    
    # move outputs into separate folders for clarity
    OUTDIR=os.path.dirname(OUTPATH)
    #moveAllFiles(outdir+"/*_SAprimers.csv", os.path.join(newoutdir, "Checkpoint_Primers"))
    #moveAllFiles(outdir+"/*_SAdimers.csv", os.path.join(newoutdir, "Checkpoint_Primers"))
    moveAllFiles(OUTDIR+"/*_primers.csv", os.path.join(OUTDIR, "Final_Primers"))
    moveAllFiles(OUTDIR+"/*_dimers.csv", os.path.join(OUTDIR, "Final_Dimers"))
    moveAllFiles(OUTDIR+"/*_ASA_costs.csv", os.path.join(OUTDIR, "Trace_Dimer_Load"))
    moveAllFiles(OUTDIR+"/*_DimerLoad.png", os.path.join(OUTDIR, "Plots_Dimer_Load"))



def moveAllFiles(filegrep, dest):
    if os.path.exists(dest):
        shutil.rmtree(dest)
    if not os.path.exists(dest):
        os.mkdir(dest)
    filelist = glob.glob(filegrep)
    for f in filelist:
        # try/except prevents errors and overwriting if same destname already exists
        try:
            shutil.move(f, dest)
        except Exception:
            print("Failed to move "+f)
            pass


class OptimizationWarning(Exception):
    pass
        

### Set up timeout exception behavior
class TimeoutException(Exception):   # Custom exception class
    pass

def timeout_handler(signum, frame):   # Custom signal handler
    raise TimeoutException


# Change the behavior of SIGALRM
signal.signal(signal.SIGALRM, timeout_handler)