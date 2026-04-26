#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 21:01:43 2026

@author: maggiehallerud
"""

# load modules
import os
import csv
import sys
import string
import primer3
import pandas as pd
import argparse

# import additional sub-modules
sys.path.append(os.path.dirname(__file__))
from CSVtoFasta import main as CSV2FASTA
from add_keeplist_to_fasta import main as AddKeeplist2FASTA



def main(TEMPLATES, OUTPATH, Tm_LIMIT=45, dG_HAIRPINS=-2, dG_END_LIMIT=-4, 
         dG_MID_LIMIT=-8, KEEPLIST=None, ENABLE_BROAD=False, SETTINGS=None):
    """
    TEMPLATES : CSV
    OUTPATH : Filename for output CSV
    Tm_LIMIT : Melting temperature threshold limit for dimers.
    dG_HAIRPINS : Delta G threshold for hairpin structures. (Default: -2000)
    dG_END_LIMIT : Delta G threshold for 3' end dimers. (Default: -4000)
    dG_MID_LIMIT : Delta G threshold for other dimers. (Default: -8000)
    KEEPLIST : Fasta file of primer pairs that must be included in final panel.
    ENABLE_BROAD : Use broader settings if no primers designed? (Default: False)
    SETTINGS : Primer3 Settings in dictionary format i.e. {'SETTINGNAME': 123, 'SETTINGNAME': "string",...}
    ------
    Output: CSV of filtered primer sequences and details
    """
    
    # read in templates as pandas dictionary
    templates = pd.read_csv(TEMPLATES)
    
    # Raise error if names aren't SEQUENCE_ID, SEQUENCE_TEMPLATE, SEQUENCE_TARGET
    if not all([x in templates.keys() for x in ['SEQUENCE_ID', 'SEQUENCE_TEMPLATE', 'SEQUENCE_TARGET']]):
        raise InputError("TEMPLATES must include fields named SEQUENCE_ID, SEQUENCE_TEMPLATE, & SEQUENCE TARGET.")
    
    # Raise error if SEQUENCE_IDs aren't unique
    if len(set(templates['SEQUENCE_ID'])) != len(templates):
        counts = templates['SEQUENCE_ID'].value_counts()
        duplicates = counts.loc[counts>1]
        print(duplicates)
        raise InputError("TEMPLATES INCLUDE NON-UNIQUE SEQUENCEID VALUES!"+
                         "   (Non-unique values and counts printed above) ")
    
    # convert to dictionary
    td = templates.set_index("SEQUENCE_ID").T.to_dict()
    
    # Check KEEPLIST_FA input before proceeding
    if KEEPLIST is not None:
        if not os.path.exists(KEEPLIST):
            raise InputError("KEEPLIST not found!")
        #keeplist = 
    
    # convert deltaG thresholds from kcal/mol to cal/mol
    dG_HAIRPINS = dG_HAIRPINS*1000
    dG_END_LIMIT = dG_END_LIMIT*1000
    dG_MID_LIMIT = dG_MID_LIMIT*1000
    
    # set up primer3 settings
    PRIMER3_SETTINGS = {'SEQUENCE_OVERHANG_LEFT': "tcgtcggcagcgtcagatgtgtataagagacag",
    'SEQUENCE_OVERHANG_RIGHT': "gtctcgtgggctcggagatgtgtataagagacag",
    'PRIMER_ANNEALING_TEMP': 52.0,
    'PRIMER_DMSO_CONC': 0.0,
    'PRIMER_DMSO_FACTOR': 0.6,
    'PRIMER_DNA_CONC': 50.0,
    'PRIMER_DNTP_CONC': 0.25,
    'PRIMER_EXPLAIN_FLAG': 1,
    'PRIMER_FIRST_BASE_INDEX': 1,
    'PRIMER_FORMAMIDE_CONC': 0.0,
    'PRIMER_GC_CLAMP': 1,
    'PRIMER_INSIDE_PENALTY': -1.0,
    'PRIMER_LIBERAL_BASE': 1,
    'PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS': 0,
    'PRIMER_LOWERCASE_MASKING': 0,
    'PRIMER_MAX_BOUND': 110.0,
    'PRIMER_MAX_END_GC': 4,
    'PRIMER_MAX_END_STABILITY': 6.0,
    'PRIMER_MAX_GC': 70,
    'PRIMER_MAX_HAIRPIN_TH': 1000.0,
    'PRIMER_MAX_LIBRARY_MISPRIMING': 12.0,
    'PRIMER_MAX_NS_ACCEPTED': 0,
    'PRIMER_MAX_POLY_X': 4,
    'PRIMER_MAX_SELF_ANY': 1000.0,
    'PRIMER_MAX_SELF_ANY_TH': 1000.0,
    'PRIMER_MAX_SELF_END': 1000.0,
    'PRIMER_MAX_SELF_END_TH': 1000.0,
    'PRIMER_MAX_SIZE': 26,
    'PRIMER_MAX_TEMPLATE_MISPRIMING': 12.0,
    'PRIMER_MAX_TEMPLATE_MISPRIMING_TH': 47.0,
    'PRIMER_MAX_TM': 63.0,
    'PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION': 4,
    'PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION': 7,
    'PRIMER_MIN_BOUND': -10.0,
    'PRIMER_MIN_END_QUALITY': 0,
    'PRIMER_MIN_GC': 30,
    'PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE': 3,
    'PRIMER_MIN_QUALITY': 0,
    'PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE': 3,
    'PRIMER_MIN_SIZE': 18,
    'PRIMER_MIN_TM': 57.0,
    'PRIMER_MISPRIMING_LIBRARY': "",
    'PRIMER_NUM_RETURN': 20,
    'PRIMER_OPT_BOUND': 97.0,
    'PRIMER_OPT_GC_PERCENT': 50.0,
    'PRIMER_OPT_SIZE': 20,
    'PRIMER_OPT_TM': 60.0,
    'PRIMER_OUTSIDE_PENALTY': 0.0,
    'PRIMER_PAIR_MAX_COMPL_ANY': 1000.0,
    'PRIMER_PAIR_MAX_COMPL_ANY_TH': 1000.0,
    'PRIMER_PAIR_MAX_COMPL_END': 1000.0,
    'PRIMER_PAIR_MAX_COMPL_END_TH': 1000.0,
    'PRIMER_PAIR_MAX_DIFF_TM': 3.0,
    'PRIMER_PAIR_MAX_LIBRARY_MISPRIMING': 24.0,
    'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING': 24.0,
    'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH': 47.0,
    'PRIMER_PAIR_WT_COMPL_ANY': 0.0,
    'PRIMER_PAIR_WT_COMPL_ANY_TH': 0.0,
    'PRIMER_PAIR_WT_COMPL_END': 0.0,
    'PRIMER_PAIR_WT_COMPL_END_TH': 0.0,
    'PRIMER_PAIR_WT_DIFF_TM': 0.0,
    'PRIMER_PAIR_WT_IO_PENALTY': 0.0,
    'PRIMER_PAIR_WT_LIBRARY_MISPRIMING': 0.0,
    'PRIMER_PAIR_WT_PRODUCT_SIZE_GT': 0.0,
    'PRIMER_PAIR_WT_PRODUCT_SIZE_LT': 0.0,
    'PRIMER_PAIR_WT_PRODUCT_TM_GT': 0.0,
    'PRIMER_PAIR_WT_PRODUCT_TM_LT': 0.0,
    'PRIMER_PAIR_WT_PR_PENALTY': 1.0,
    'PRIMER_PAIR_WT_TEMPLATE_MISPRIMING': 0.0,
    'PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH': 0.0,
    'PRIMER_PICK_ANYWAY': 0,
    'PRIMER_PICK_INTERNAL_OLIGO': 0,
    'PRIMER_PICK_LEFT_PRIMER': 1,
    'PRIMER_PICK_RIGHT_PRIMER': 1,
    'PRIMER_PRODUCT_MAX_TM': 1000000.0,
    'PRIMER_PRODUCT_MIN_TM': -1000000.0,
    'PRIMER_PRODUCT_OPT_SIZE': 100,
    'PRIMER_PRODUCT_OPT_TM': 0.0,
    'PRIMER_PRODUCT_SIZE_RANGE': "70-120",
    'PRIMER_QUALITY_RANGE_MAX': 100,
    'PRIMER_QUALITY_RANGE_MIN': 0,
    'PRIMER_SALT_CORRECTIONS': 1,
    'PRIMER_SALT_DIVALENT': 3.8,
    'PRIMER_SALT_MONOVALENT': 50.0,
    'PRIMER_SECONDARY_STRUCTURE_ALIGNMENT': 1,
    'PRIMER_SEQUENCING_ACCURACY': 20,
    'PRIMER_SEQUENCING_INTERVAL': 250,
    'PRIMER_SEQUENCING_LEAD': 50,
    'PRIMER_SEQUENCING_SPACING': 500,
    'PRIMER_TASK': "generic",
    'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,
    'PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT': 0,
    'PRIMER_TM_FORMULA': 1,
    'PRIMER_WT_BOUND_GT': 0.0,
    'PRIMER_WT_BOUND_LT': 0.0,
    'PRIMER_WT_END_QUAL': 0.0,
    'PRIMER_WT_END_STABILITY': 0.0,
    'PRIMER_WT_GC_PERCENT_GT': 0.0,
    'PRIMER_WT_GC_PERCENT_LT': 0.0,
    'PRIMER_WT_HAIRPIN_TH': 0.0,
    'PRIMER_WT_LIBRARY_MISPRIMING': 0.0,
    'PRIMER_WT_NUM_NS': 0.0,
    'PRIMER_WT_POS_PENALTY': 0.0,
    'PRIMER_WT_SELF_ANY': 0.0,
    'PRIMER_WT_SELF_ANY_TH': 0.0,
    'PRIMER_WT_SELF_END': 0.0,
    'PRIMER_WT_SELF_END_TH': 0.0,
    'PRIMER_WT_SEQ_QUAL': 0.0,
    'PRIMER_WT_SIZE_GT': 1.0,
    'PRIMER_WT_SIZE_LT': 1.0,
    'PRIMER_WT_TEMPLATE_MISPRIMING': 0.0,
    'PRIMER_WT_TEMPLATE_MISPRIMING_TH': 0.0,
    'PRIMER_WT_TM_GT': 1.0,
    'PRIMER_WT_TM_LT': 1.0,
    'P3_PILE_FLAG': 1,
    'NONSENSE': 'ABCD'}
    
    BROAD_SETTINGS = {'SEQUENCE_OVERHANG_LEFT': "tcgtcggcagcgtcagatgtgtataagagacag",
    'SEQUENCE_OVERHANG_RIGHT': "gtctcgtgggctcggagatgtgtataagagacag",
    'PRIMER_ANNEALING_TEMP': 52.0,
    'PRIMER_DMSO_CONC': 0.0,
    'PRIMER_DMSO_FACTOR': 0.6,
    'PRIMER_DNA_CONC': 50.0,
    'PRIMER_DNTP_CONC': 0.25,
    'PRIMER_EXPLAIN_FLAG': 1,
    'PRIMER_FIRST_BASE_INDEX': 1,
    'PRIMER_FORMAMIDE_CONC': 0.0,
    'PRIMER_GC_CLAMP': 0,
    'PRIMER_INSIDE_PENALTY': -1.0,
    'PRIMER_LIBERAL_BASE': 1,
    'PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS': 0,
    'PRIMER_LOWERCASE_MASKING': 0,
    'PRIMER_MAX_BOUND': 110.0,
    'PRIMER_MAX_END_GC': 5,
    'PRIMER_MAX_END_STABILITY': 6.0,
    'PRIMER_MAX_GC': 80,
    'PRIMER_MAX_HAIRPIN_TH': 1000.00,
    'PRIMER_MAX_LIBRARY_MISPRIMING': 1000.00,
    'PRIMER_MAX_NS_ACCEPTED': 0,
    'PRIMER_MAX_POLY_X': 5,
    'PRIMER_MAX_SELF_ANY': 1000.00,
    'PRIMER_MAX_SELF_ANY_TH': 1000.00,
    'PRIMER_MAX_SELF_END': 1000.00,
    'PRIMER_MAX_SELF_END_TH': 1000.00,
    'PRIMER_MAX_SIZE': 26,
    'PRIMER_MAX_TEMPLATE_MISPRIMING': 12.00,
    'PRIMER_MAX_TEMPLATE_MISPRIMING_TH': 47.00,
    'PRIMER_MAX_TM': 63.0,
    'PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION': 4,
    'PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION': 7,
    'PRIMER_MIN_BOUND': -10.0,
    'PRIMER_MIN_END_QUALITY': 0,
    'PRIMER_MIN_GC': 20,
    'PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE': 3,
    'PRIMER_MIN_QUALITY': 0,
    'PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE': 3,
    'PRIMER_MIN_SIZE': 18,
    'PRIMER_MIN_TM': 57.0,
    'PRIMER_MISPRIMING_LIBRARY': "",
    'PRIMER_NUM_RETURN': 20,
    'PRIMER_OPT_BOUND': 97.0,
    'PRIMER_OPT_GC_PERCENT': 50.0,
    'PRIMER_OPT_SIZE': 20,
    'PRIMER_OPT_TM': 60.0,
    'PRIMER_OUTSIDE_PENALTY': 0.0,
    'PRIMER_PAIR_MAX_COMPL_ANY': 1000.00,
    'PRIMER_PAIR_MAX_COMPL_ANY_TH': 1000.00,
    'PRIMER_PAIR_MAX_COMPL_END': 1000.00,
    'PRIMER_PAIR_MAX_COMPL_END_TH': 1000.00,
    'PRIMER_PAIR_MAX_DIFF_TM': 3.0,
    'PRIMER_PAIR_MAX_LIBRARY_MISPRIMING': 24.00,
    'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING': 24.00,
    'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH': 47.00,
    'PRIMER_PAIR_WT_COMPL_ANY': 0.0,
    'PRIMER_PAIR_WT_COMPL_ANY_TH': 0.0,
    'PRIMER_PAIR_WT_COMPL_END': 0.0,
    'PRIMER_PAIR_WT_COMPL_END_TH': 0.0,
    'PRIMER_PAIR_WT_DIFF_TM': 0.0,
    'PRIMER_PAIR_WT_IO_PENALTY': 0.0,
    'PRIMER_PAIR_WT_LIBRARY_MISPRIMING': 0.0,
    'PRIMER_PAIR_WT_PRODUCT_SIZE_GT': 0.0,
    'PRIMER_PAIR_WT_PRODUCT_SIZE_LT': 0.0,
    'PRIMER_PAIR_WT_PRODUCT_TM_GT': 0.0,
    'PRIMER_PAIR_WT_PRODUCT_TM_LT': 0.0,
    'PRIMER_PAIR_WT_PR_PENALTY': 1.0,
    'PRIMER_PAIR_WT_TEMPLATE_MISPRIMING': 0.0,
    'PRIMER_PAIR_WT_TEMPLATE_MISPRING_TH': 0.0,
    'PRIMER_PICK_ANYWAY': 0,
    'PRIMER_PICK_INTERNAL_OLIGO': 0,
    'PRIMER_PICK_LEFT_PRIMER': 1,
    'PRIMER_PICK_RIGHT_PRIMER': 1,
    'PRIMER_PRODUCT_MAX_TM': 1000000.0,
    'PRIMER_PRODUCT_MIN_TM': -1000000.0,
    'PRIMER_PRODUCT_OPT_SIZE': 100,
    'PRIMER_PRODUCT_OPT_TM': 0.0,
    'PRIMER_PRODUCT_SIZE_RANGE': "70-120",
    'PRIMER_QUALITY_RANGE_MAX': 100,
    'PRIMER_QUALITY_RANGE_MIN': 0,
    'PRIMER_SALT_CORRECTIONS': 1,
    'PRIMER_SALT_DIVALENT': 3.8,
    'PRIMER_SALT_MONOVALENT': 50.0,
    'PRIMER_SELF_ANY': 8.00,
    'PRIMER_SELF_END': 3.00,
    'PRIMER_SECONDARY_STRUCTURE_ALIGNMENT': 1,
    'PRIMER_SEQUENCING_ACCURACY': 20,
    'PRIMER_SEQUENCING_INTERVAL': 250,
    'PRIMER_SEQUENCING_LEAD': 50,
    'PRIMER_SEQUENCING_SPACING': 500,
    'PRIMER_TASK': "generic",
    'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,
    'PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT': 0,
    'PRIMER_TM_FORMULA': 1,
    'PRIMER_WT_BOUND_GT': 0.0,
    'PRIMER_WT_BOUND_LT': 0.0,
    'PRIMER_WT_END_QUAL': 0.0,
    'PRIMER_WT_END_STABILITY': 0.0,
    'PRIMER_WT_GC_PERCENT_GT': 0.0,
    'PRIMER_WT_GC_PERCENT_LT': 0.0,
    'PRIMER_WT_HAIRPIN_TH': 0.0,
    'PRIMER_WT_LIBRARY_MISPRIMING': 0.0,
    'PRIMER_WT_NUM_NS': 0.0,
    'PRIMER_WT_POS_PENALTY': 0.0,
    'PRIMER_WT_SELF_ANY': 0.0,
    'PRIMER_WT_SELF_ANY_TH': 0.0,
    'PRIMER_WT_SELF_END': 0.0,
    'PRIMER_WT_SELF_END_TH': 0.0,
    'PRIMER_WT_SEQ_QUAL': 0.0,
    'PRIMER_WT_SIZE_GT': 1.0,
    'PRIMER_WT_SIZE_LT': 1.0,
    'PRIMER_WT_TEMPLATE_MISPRIMING': 0.0,
    'PRIMER_WT_TEMPLATE_MISPRIMING_TH': 0.0,
    'PRIMER_WT_TM_GT': 1.0,
    'PRIMER_WT_TM_LT': 1.0,
    'P3_PILE_FLAG': 1 }
    
    if SETTINGS is not None:
        for k, v in SETTINGS.items():
            if k in PRIMER3_SETTINGS:
                # check type of setting
                if type(PRIMER3_SETTINGS[k]) != type(SETTINGS[k]):
                    print("Primer3 setting "+str(k)+" could not be adjusted- Format should be "+str(type(PRIMER3_SETTINGS[k])))
                # update otherwise.
                else:                    
                    try:
                        PRIMER3_SETTINGS[k] = v
                        BROAD_SETTINGS[k] = v
                    except Exception:
                        print("Primer3 setting "+str(k)+" could not be adjusted. "+
                                 "Default value ("+str(PRIMER3_SETTINGS[k])+") used.")
            else:
                print("Primer3 setting "+str(k)+" not found- are you sure this is a primer3 setting?")    
                try:
                    PRIMER3_SETTINGS[k] = v
                    BROAD_SETTINGS[k] = v
                except Exception:
                    pass
    
    Ladapt = PRIMER3_SETTINGS['SEQUENCE_OVERHANG_LEFT'].upper()
    Radapt = PRIMER3_SETTINGS['SEQUENCE_OVERHANG_RIGHT'].upper()
    
    # set up empty tuple to hold filtered output
    filtered_primers = []
    filtered_primers.append(["PrimerID","LocusID","PrimerPair","Direction","Sequence",
                             "AmpliconSize","StartBP","Length","AnnealingTempC",
                             "GC%","PropBound","EndStability","Penalty","TemplMispriming"])
    
    # set up arrays for progress tracking
    passed_ids = []
    failed_ids = []
    tot_pairs = 0

    # set up translation table for replacing all punctuation in input
    translator = str.maketrans('', '', string.punctuation.replace("_","").replace(":",""))
    
    # loop through each template, design & filter primers
    row=0
    for k in td.keys():
        # progress tracking
        row+=1
        if row%100 == 0:
            print("      primers designed for "+str(row)+" sequences")
        
        # extract sequence ID, template, target info
        seqid = str(k).translate(translator)
        # design primers for template
        try:
            out = primer3.bindings.design_primers(seq_args={
                        'SEQUENCE_ID': k,
                        'SEQUENCE_TEMPLATE': td[k]['SEQUENCE_TEMPLATE'],
                        'SEQUENCE_TARGET': [int(x) for x in td[k]['SEQUENCE_TARGET'].split(",")]}, 
                    global_args=PRIMER3_SETTINGS)
        except Exception as err:
            print("PRIMER DESIGN FAILED FOR "+k+" WITH THE FOLLOWING ERROR:")
            print(err)
            

        # save full primer3 output to file
        # remove lists- this info is all redundant
        #out.pop('PRIMER_PAIR')
        #out.pop('PRIMER_LEFT')
        #out.pop('PRIMER_RIGHT')
        #with open() as file:
        #    for k, v in out.items():
        #        file.write(f"{k}={v}\n")
        
        # try running with broader settings if none found
        pairs = out['PRIMER_PAIR_NUM_RETURNED']
        if pairs == 0:
            if ENABLE_BROAD:
                try:
                    out = primer3.bindings.design_primers(seq_args={
                                'SEQUENCE_ID': k,
                                'SEQUENCE_TEMPLATE': td[k]['SEQUENCE_TEMPLATE'],
                                'SEQUENCE_TARGET': [int(x) for x in td[k]['SEQUENCE_TARGET'].split(",")]}, 
                            global_args=BROAD_SETTINGS)
                except Exception as err:
                    print("PRIMER DESIGN FAILED FOR "+k+" WITH THE FOLLOWING ERROR:")
                    print(err)
        
        # loop through primer pairs & filter dimers     
        passed = 0 # track number pairs passed per template
        pairs = out['PRIMER_PAIR_NUM_RETURNED']
        if pairs > 0:
            for N in range(pairs):
                tot_pairs +=pairs # count overall # primer pairs designed
                N = str(N)
                # pull primer sequences for this pair
                FWseq = out['PRIMER_LEFT_'+N+'_SEQUENCE']
                RVseq = out['PRIMER_RIGHT_'+N+'_SEQUENCE']
                # calculate hairpins - primer + adapter
                hairpinL = primer3.bindings.calc_hairpin(FWseq, 
                                                         mv_conc=PRIMER3_SETTINGS['PRIMER_SALT_MONOVALENT'],
                                                         dv_conc=PRIMER3_SETTINGS['PRIMER_SALT_DIVALENT'], 
                                                         dntp_conc=PRIMER3_SETTINGS['PRIMER_DNTP_CONC'],
                                                         dna_conc=PRIMER3_SETTINGS['PRIMER_DNA_CONC'])
                hairpinR = primer3.bindings.calc_hairpin(RVseq, 
                                                         mv_conc=PRIMER3_SETTINGS['PRIMER_SALT_MONOVALENT'],
                                                         dv_conc=PRIMER3_SETTINGS['PRIMER_SALT_DIVALENT'], 
                                                         dntp_conc=PRIMER3_SETTINGS['PRIMER_DNTP_CONC'],
                                                         dna_conc=PRIMER3_SETTINGS['PRIMER_DNA_CONC'])
                # calculate hairpins - primer only
                hairpinPL = primer3.bindings.calc_hairpin(FWseq.replace(Ladapt, ''), 
                                                         mv_conc=PRIMER3_SETTINGS['PRIMER_SALT_MONOVALENT'],
                                                         dv_conc=PRIMER3_SETTINGS['PRIMER_SALT_DIVALENT'], 
                                                         dntp_conc=PRIMER3_SETTINGS['PRIMER_DNTP_CONC'],
                                                         dna_conc=PRIMER3_SETTINGS['PRIMER_DNA_CONC'])
                hairpinPR = primer3.bindings.calc_hairpin(RVseq.replace(Radapt, ''), 
                                                         mv_conc=PRIMER3_SETTINGS['PRIMER_SALT_MONOVALENT'],
                                                         dv_conc=PRIMER3_SETTINGS['PRIMER_SALT_DIVALENT'], 
                                                         dntp_conc=PRIMER3_SETTINGS['PRIMER_DNTP_CONC'],
                                                         dna_conc=PRIMER3_SETTINGS['PRIMER_DNA_CONC'])
                # calculate dimers- most likely
                homodimerL = primer3.bindings.calc_homodimer(FWseq, 
                                                         mv_conc=PRIMER3_SETTINGS['PRIMER_SALT_MONOVALENT'],
                                                         dv_conc=PRIMER3_SETTINGS['PRIMER_SALT_DIVALENT'], 
                                                         dntp_conc=PRIMER3_SETTINGS['PRIMER_DNTP_CONC'],
                                                         dna_conc=PRIMER3_SETTINGS['PRIMER_DNA_CONC'])
                homodimerR = primer3.bindings.calc_homodimer(RVseq, 
                                                         mv_conc=PRIMER3_SETTINGS['PRIMER_SALT_MONOVALENT'],
                                                         dv_conc=PRIMER3_SETTINGS['PRIMER_SALT_DIVALENT'], 
                                                         dntp_conc=PRIMER3_SETTINGS['PRIMER_DNTP_CONC'],
                                                         dna_conc=PRIMER3_SETTINGS['PRIMER_DNA_CONC'])
                dimers = primer3.bindings.calc_heterodimer(FWseq, RVseq, 
                                                         mv_conc=PRIMER3_SETTINGS['PRIMER_SALT_MONOVALENT'],
                                                         dv_conc=PRIMER3_SETTINGS['PRIMER_SALT_DIVALENT'], 
                                                         dntp_conc=PRIMER3_SETTINGS['PRIMER_DNTP_CONC'],
                                                         dna_conc=PRIMER3_SETTINGS['PRIMER_DNA_CONC'])
                # calculate end dimers
                enddimerL = primer3.bindings.calc_end_stability(FWseq, FWseq, 
                                                         mv_conc=PRIMER3_SETTINGS['PRIMER_SALT_MONOVALENT'],
                                                         dv_conc=PRIMER3_SETTINGS['PRIMER_SALT_DIVALENT'], 
                                                         dntp_conc=PRIMER3_SETTINGS['PRIMER_DNTP_CONC'],
                                                         dna_conc=PRIMER3_SETTINGS['PRIMER_DNA_CONC'])
                enddimerR = primer3.bindings.calc_end_stability(RVseq, RVseq, 
                                                         mv_conc=PRIMER3_SETTINGS['PRIMER_SALT_MONOVALENT'],
                                                         dv_conc=PRIMER3_SETTINGS['PRIMER_SALT_DIVALENT'], 
                                                         dntp_conc=PRIMER3_SETTINGS['PRIMER_DNTP_CONC'],
                                                         dna_conc=PRIMER3_SETTINGS['PRIMER_DNA_CONC'])
                enddimerLR = primer3.bindings.calc_end_stability(FWseq, RVseq, 
                                                         mv_conc=PRIMER3_SETTINGS['PRIMER_SALT_MONOVALENT'],
                                                         dv_conc=PRIMER3_SETTINGS['PRIMER_SALT_DIVALENT'], 
                                                         dntp_conc=PRIMER3_SETTINGS['PRIMER_DNTP_CONC'],
                                                         dna_conc=PRIMER3_SETTINGS['PRIMER_DNA_CONC'])
                # run tests for each structure
                tests = 0 #set counter
                tests = testStructure(hairpinL, Tm_LIMIT, dG_HAIRPINS, tests)
                tests = testStructure(hairpinR, Tm_LIMIT, dG_HAIRPINS, tests)
                tests = testStructure(hairpinPL, Tm_LIMIT, dG_HAIRPINS, tests)
                tests = testStructure(hairpinPR, Tm_LIMIT, dG_HAIRPINS, tests)
                tests = testStructure(homodimerL, Tm_LIMIT, dG_MID_LIMIT, tests)
                tests = testStructure(homodimerR, Tm_LIMIT, dG_MID_LIMIT, tests)
                tests = testStructure(dimers, Tm_LIMIT, dG_MID_LIMIT, tests)
                tests = testStructure(enddimerL, Tm_LIMIT, dG_END_LIMIT, tests)
                tests = testStructure(enddimerR, Tm_LIMIT, dG_END_LIMIT, tests)
                tests = testStructure(enddimerLR, Tm_LIMIT, dG_END_LIMIT, tests)
                # export info if all tests were passed
                if tests == 10:
                    ampSize = out['PRIMER_PAIR_'+N+'_PRODUCT_SIZE']
                    for DIR in ['LEFT','RIGHT']:
                        # extract useful info
                        pos = out['PRIMER_'+DIR+'_'+N]
                        tm = out['PRIMER_'+DIR+'_'+N+'_TM']
                        bound = out['PRIMER_'+DIR+'_'+N+'_BOUND']
                        gc = out['PRIMER_'+DIR+'_'+N+'_GC_PERCENT']
                        endstab = out['PRIMER_'+DIR+'_'+N+'_END_STABILITY']
                        penalty =out['PRIMER_'+DIR+'_'+N+'_PENALTY']
                        misprime = out['PRIMER_'+DIR+'_'+N+'_TEMPLATE_MISPRIMING']
                        if DIR=="LEFT":
                            DIR = "FWD"
                            seq=FWseq
                        if DIR=="RIGHT":
                            DIR = "REV"
                            seq = RVseq
                        # add to filtered primers list
                        filtered_primers.append([seqid+"."+N+"."+DIR, # primer ID
                                                 seqid, # locus ID
                                                 N, # primer pair #
                                                 DIR, # primer direction
                                                 seq.upper(), # primer sequence (capitalized)
                                                 ampSize, #amplicon size
                                                 pos[0], # primer position (start BP)
                                                 pos[1], # primer position (length)
                                                 round(tm,1), # primer annealing temp
                                                 round(gc,1), # gc content
                                                 round(bound), # primer proportion bound
                                                 round(endstab,2), # end stability
                                                 round(penalty,2), #primer3 penalty
                                                 round(misprime,1)])#misprime rate
                    # progress tracking
                    passed +=1 #keep track of # pairs passing per template
            # track whether template had passing primer pairs or not
            if passed > 0:
                passed_ids.append(seqid)
            else:
                failed_ids.append(seqid)
    
    # Print number loci in filtered output
    print("# Input loci: "+str(len(templates)))
    print("# Loci with primers passing filtering: "+str(len(passed_ids)))
    print("# Primer pairs designed: "+str(tot_pairs))
    print("# Primer pairs passing filtering: "+str(len(filtered_primers)))
    if len(failed_ids)>0:
        print("Filtering failed for the following loci:")
        for l in failed_ids:
            print("          "+l)

    # Export filtered primers as CSV
    if len(filtered_primers)>1:
        OUTCSV= OUTPATH + '.csv'
        with open(OUTCSV, 'w', newline="\n") as file:
            writer = csv.writer(file)
            for row in filtered_primers:
                writer.writerow(row)
    
        # Export filtered primers as FASTA file
        CSV2FASTA(OUTCSV, OUTPATH+".fa")
        
        # if keeplist provided, add 
        if KEEPLIST is not None:
            try: 
                AddKeeplist2FASTA(OUTPATH+".fa", KEEPLIST)
            except Exception:
                print("KEEPLIST could not be added to filtered primers FASTA")                

    # raise Exception here if no primer passed filtering since pointless to continue
    else:
        raise OutputError("NO PRIMER PAIRS PASSED FILTERING- STOPPING AT DESIGN/FILTER STEP!")
    


def testStructure(structure, TM_LIMIT, DG_LIMIT, TESTS):
    dg = structure.dg
    tm = structure.tm
    if dg>DG_LIMIT or tm<TM_LIMIT:
        TESTS+=1
    return(TESTS)



class InputError(Exception):
    pass



class OutputError(Exception):
    pass



def parse_args():
    # initialize argparser
    parser = argparse.ArgumentParser()
    # add required arguments
    parser.add_argument("-t", "--templates", type=str, required=True)
    parser.add_argument("-o", "--out", type=str, required=True)
    # add optional arguments
    parser.add_argument("-l", "--tm_limit", type=float, default=45.0)
    parser.add_argument("-h", "--hairpins_dg", type=float, default=-2)
    parser.add_argument("-m", "--middimers_dg", type=float, default=-8)
    parser.add_argument("-e", "--enddimers_dg", type=float, default=-4)
    parser.add_argument("-k", "--keeplist", type=str, default=None)
    parser.add_argument("-b", "--enable_broad", action="store_true")
    parser.add_argument("-s", "--settings", type=dict, default=None)
    
    return parser.parse_args()



if __name__=="__main__":
    # parse command line arguments
    args = parse_args()
    # run main
    main(TEMPLATES=args.templates,
         OUTPATH=args.out, 
         Tm_LIMIT=args.tm_limit, 
         dG_HAIRPINS=args.hairpins_dg, 
         dG_END_LIMIT=args.enddimers_dg, 
         dG_MID_LIMIT=args.middimers_dg, 
         KEEPLIST=args.keeplist, 
         ENABLE_BROAD=args.enable_broad, 
         SETTINGS=args.settings)
