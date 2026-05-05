#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TITLE: PrimerBLAST Thermodynamic Calculations
Purpose: Calculates Tm and delta G of off-target primer-binding sites
        predicted by PRIMER-BLAST

Created on Mon May  4 13:35:50 2026
@author: maggiehallerud
"""
# load dependencies
import os
import primer3
import argparse
import pandas as pd
import warnings
pd.options.mode.chained_assignment = None #default 'warn' (none avoids unnecessary warnings)



def main(INFILE, OUTFILE, ANNEAL_TEMP=52.0, MV_CONC=50, DV_CONC=3.8, DNTP_CONC=0.25, DNA_CONC=50):
    """
    INFILE : PRIMER-BLAST outputs.
    OUTFILE : Path to save CSV.
    ANNEAL_TEMP : PCR annealing temp (Celsius). Default: 52
    MV_CONC : Monovalent salt concentration in PCR. Default: 50
    DV_CONC : Divalent salt concentration in PCR. Default: 3.8
    DNTP_CONC : Primer concentration in PCR. Default: 0.25
    DNA_CONC : DNA concentration in PCR. Default: 50
    ----------
    """
    # check inputs
    if not os.path.exists(INFILE):
        raise InputError("INFILE file not found!")
    
    # read in primer-blast inputs
    inputs = pd.read_csv(INFILE)
    
    # re-calculate primer start/stop locations based on product position
    inputs.forward_start = inputs['forward_start'].sub(inputs.product_start)
    inputs.forward_stop = inputs['forward_stop'].sub(inputs.product_start)
    inputs.reverse_start = inputs['reverse_start'].sub(inputs.product_start)
    inputs.reverse_stop = inputs['reverse_stop'].sub(inputs.product_start)

    # extract primer-binding sites for each off-target amplicon
    inputs['forward_binding'] = pd.Series(dtype="str")
    inputs['reverse_binding'] = pd.Series(dtype="str")
    inputs['Tm_FWD'] = pd.Series(dtype="float")
    inputs['dG_FWD'] = pd.Series(dtype="float")
    inputs['Tm_FWD_END'] = pd.Series(dtype="float")
    inputs['dG_FWD_END'] = pd.Series(dtype="float")
    inputs['Tm_REV'] = pd.Series(dtype="float")
    inputs['dG_REV'] = pd.Series(dtype="float")
    inputs['Tm_REV_END'] = pd.Series(dtype="float")
    inputs['dG_REV_END'] = pd.Series(dtype="float")
    for row in range(len(inputs)):
        fwd1 = inputs.forward_start[row]
        fwd2 = inputs.forward_stop[row]
        rev1 = inputs.reverse_start[row]
        rev2 = inputs.reverse_stop[row]
        if inputs.Sequence[row] is not None:
            # depends on whether fwd is 5' or 3' end of sequence...
            if fwd1 < fwd2:
                inputs.forward_binding[row] = inputs.loc[row,'Sequence'][fwd1:fwd2]
                inputs.reverse_binding[row] = inputs.loc[row,'Sequence'][rev2:rev1]
                FWD = inputs.FWDseq[row]
                REV = rc(inputs.REVseq[row])
            else:
                inputs.forward_binding[row] = inputs.loc[row,'Sequence'][fwd2:fwd1]
                inputs.reverse_binding[row] = inputs.loc[row,'Sequence'][rev1:rev2]
                FWD = rc(inputs.FWDseq[row])
                REV = inputs.REVseq[row]
            
            # calculate dG / Tm for each primer binding site
            fwd_binding = primer3.calc_heterodimer(inputs.forward_binding[row], FWD,
                                                   temp_c=ANNEAL_TEMP, mv_conc=MV_CONC, dv_conc=DV_CONC,
                                                   dna_conc=DNA_CONC, dntp_conc=DNTP_CONC)
            fwd_end = primer3.calc_end_stability(inputs.forward_binding[row], FWD,
                                                 temp_c=ANNEAL_TEMP, mv_conc=MV_CONC, dv_conc=DV_CONC,
                                                 dna_conc=DNA_CONC, dntp_conc=DNTP_CONC)
            rev_binding = primer3.calc_heterodimer(inputs.reverse_binding[row], REV,
                                                   temp_c=ANNEAL_TEMP, mv_conc=MV_CONC, dv_conc=DV_CONC,
                                                   dna_conc=DNA_CONC, dntp_conc=DNTP_CONC)
            rev_end = primer3.calc_end_stability(inputs.reverse_binding[row], REV,
                                                 temp_c=ANNEAL_TEMP, mv_conc=MV_CONC, dv_conc=DV_CONC,
                                                 dna_conc=DNA_CONC, dntp_conc=DNTP_CONC)
            inputs.Tm_FWD[row] = fwd_binding.tm
            inputs.dG_FWD[row] = fwd_binding.dg/1000
            inputs.Tm_REV[row] = rev_binding.tm
            inputs.dG_REV[row] = rev_binding.dg/1000
            inputs.Tm_FWD_END[row] = fwd_end.tm
            inputs.dG_FWD_END[row] = fwd_end.dg/1000
            inputs.Tm_REV_END[row] = rev_end.tm
            inputs.dG_REV_END[row] = rev_end.dg/1000

        else:
            warnings.warn("Sequence not present for row %s in INFILE", row)
                
    # save outputs
    inputs.to_csv(OUTFILE)



def rc(seq):
    complement_table = str.maketrans("ATCGatcg", "TAGCtagc")
    rc = seq.translate(complement_table)[::-1]
    return rc

    

class InputError(Exception):
    pass 



def parse_args():
    # initialize argparser
    parser = argparse.ArgumentParser()
    # add required arguments
    parser.add_argument("-i", "--infile", type=str, required=True)
    parser.add_argument("-o", "--outfile", type=str, required=True)
    # add optional arguments
    parser.add_argument("-t", "--anneal-temp", type=float, default=52)
    parser.add_argument("-m", "--mv_conc", type=float, default=50)
    parser.add_argument("-d", "--dv_conc", type=float, default=3.8)
    parser.add_argument("-p", "--dntp_conc", type=float, default=0.25)
    parser.add_argument("-c", "--dna_conc", type=float, default=50)
    return parser.parse_args()



if __name__=="__main__":
    args = parse_args()
    main(INFILE=args.infile,
         OUTFILE=args.outfile,
         ANNEAL_TEMP=args.anneal_temp,
         MV_CONC=args.mv_conc, 
         DV_CONC=args.dv_conc,
         DNTP_CONC=args.dntp_conc,
         DNA_CONC=args.dna_conc)