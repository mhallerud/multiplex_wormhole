#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 14:31:28 2024

@author: maggiehallerud
"""

import traceback
import os
import sys
import glob
import shutil
#import time
import argparse

sys.path.append(os.path.dirname(__file__))
from optimize_primers import main as optimizeMultiplex



def main(N_RUNS, PRIMER_FA, DIMER_SUMS, DIMER_TABLE, OUTPATH, N_LOCI, deltaG=False,
                          KEEPLIST=None, TIMEOUT=5, VERBOSE=False, SEED=None,
                          SIMPLE=5000, ITERATIONS=10000, BURNIN=100, DECAY_RATE=0.98, 
                          T_INIT=None, T_FINAL=None, PARTITIONS=1000, DIMER_ADJ=0.1, PROB_ADJ=2):
    # set up empty array to hold overall dimer load 
    loads = [['Run', 'TotalDimers']]
    
    # while instead of for loop allows runs that timeout to restart
    run = 1
    while run <= N_RUNS:
        # set timer- this helps to timeout runs with infinite loops in the optimization process
        #TIMEOUT = TIMEOUT*60 # convert to seconds
        #start=time.time()
        #while time.time()-start <= TIMEOUT:
            try:
                cost = optimizeMultiplex(PRIMER_FASTA=PRIMER_FA,
                                         DIMER_SUMS=DIMER_SUMS, 
                                         DIMER_TABLE=DIMER_TABLE,
                                         OUTPATH=OUTPATH+"_Run"+str(run).zfill(2),
                                         N_LOCI=N_LOCI, 
                                         KEEPLIST=KEEPLIST,
                                         deltaG=deltaG,
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
                run+=1
                print(" ")
            except Exception as e:
                print("AN ERROR OCCURRED IN RUN "+str(run))
                print(e)
                traceback.print_exc()
                print("---")
                run+=1
                continue
        #else:
        #    raise TimeoutException("in run "+str(run))
        #continue

    
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
        os.makedirs(dest, exist_ok=True)
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

#def timeout_handler(signum, frame):   # Custom signal handler
#    raise TimeoutException


# Change the behavior of SIGALRM
#signal.signal(signal.SIGALRM, timeout_handler)


def parse_args():
    # initialize argparser
    parser = argparse.ArgumentParser()
    # add required arguments
    parser.add_argument("-r", "--runs", type=int, required=True)
    parser.add_argument("-f", "--primer_fasta", type=str, required=True)
    parser.add_argument("-d", "--dimer_sums", type=str, required=True)
    parser.add_argument("-t", "--dimer_table", type=str, required=True)
    parser.add_argument("-o", "--outpath", type=str, required=True)
    parser.add_argument("-n", "--nloci", type=int, required=True)
    # add optional arguments
    parser.add_argument("-k", "--keeplist", type=str, default=None)
    parser.add_argument("-z", "--seed", type=str, default=None)
    parser.add_argument("-s", "--simple", type=int, default=5000)
    parser.add_argument("-i", "--iter", type=int, default=10000)
    parser.add_argument("-b", "--burnin", type=int, default=100)
    parser.add_argument("-y", "--decay_rate", type=float, default=0.95)
    parser.add_argument("-x", "--temp_init", type=float, default=None)
    parser.add_argument("-l", "--temp_final", type=float, default=0.1)
    parser.add_argument("-n", "--partitions", type=int, default=1000)
    parser.add_argument("-w", "--dimer_adj", type=float, default=0.1)
    parser.add_argument("-a", "--prob_adj", type=float, default=2)
    parser.add_argument("-u", "--timeout", type=float, default=5)
    # add flags
    parser.add_argument("-g", "--deltaG", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-m", "--makeplot", action="store_true")
    return parser.parse_args()



if __name__=="__main__":
    # parse command line arguments
    args = parse_args()
    # run main
    main(N_RUNS = args.runs, 
         PRIMER_FA = args.primer_fasta, 
         DIMER_SUMS = args.dimer_sums, 
         DIMER_TABLE = args.dimer_table, 
         OUTPATH = args.outpath, 
         N_LOCI = args.nloci, 
         deltaG = args.deltaG,
         KEEPLIST = args.keeplist, 
         TIMEOUT = args.timeout, 
         VERBOSE = args.verbose, 
         SEED = args.seed,
         SIMPLE = args.simple, 
         ITERATIONS = args.iterations, 
         BURNIN = args.burnin, 
         DECAY_RATE = args.decay_rate,
         T_INIT = args.temp_init, 
         T_FINAL = args.temp_final, 
         PARTITIONS = args.partitions, 
         DIMER_ADJ = args.dimer_adj, 
         PROB_ADJ = args.prob_adj)
