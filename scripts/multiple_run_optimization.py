#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 14:31:28 2024

@author: maggiehallerud
"""

import signal
from optimize_primers import main as optimizeMultiplex

def multipleOptimizations(N_RUNS, PRIMER_FA, DIMER_SUMS, DIMER_TABLE, OUTPATH, N_LOCI, WHITELIST=None, TIMEOUT=10):
    # set up empty array to hold overall dimer load 
    loads = [['Run', 'TotalDimers']]
    
    # while instead of for loop allows runs that timeout to restart
    run = 1
    while run <= N_RUNS:
        signal.alarm(TIMEOUT*60) # set timer- this helps to timeout runs with infinite loops in the optimization process
        try:
            if WHITELIST is None:
                cost = optimizeMultiplex(PRIMER_FA, DIMER_SUMS, DIMER_TABLE, OUTPATH+"_Run"+str(run), N_LOCI, WHITELIST=None)
            else:
                cost = optimizeMultiplex(PRIMER_FA, DIMER_SUMS, DIMER_TABLE, OUTPATH+"_Run"+str(run), N_LOCI, WHITELIST)
            loads.append([str(run), str(cost)])
        except TimeoutException:
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



### Set up timeout exception behavior
class TimeoutException(Exception):   # Custom exception class
    pass

def timeout_handler(signum, frame):   # Custom signal handler
    raise TimeoutException

# Change the behavior of SIGALRM
signal.signal(signal.SIGALRM, timeout_handler)