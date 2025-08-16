#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 11:06:54 2024

@author: maggiehallerud
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 18:36:10 2024

@author: maggiehallerud
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: TABULATE MFEprimer DIMERS
Created on Wed Aug  2 11:57:09 2023
@author: maggiehallerud

Dependencies: pandas (developed with version 1.4.4)

Purpose: Converts text output from MFEprimer dimer function into a pairwise interaction table
        and sum of primer interactions per primer and primer pair
"""

# load dependencies
import sys
import gc
import pandas as pd #version 1.4.4
import itertools #paralellized filtering, etc.
#from operator import truth # converts anything>=1 to True and =0 to False
#import multiprocessing



def main(ALL_DIMERS, END_DIMERS, OUTPATH, OUTPRIMERPATH="False"):
    """
    ALL_DIMERS : Filepath
        Output from MFEprimer dimer
    END_DIMERS : Filepath
        Output from MFEprimer dimer -p
    OUTPATH : Filepath
        Filepath prefix for all primer pair outputs
    OUTPRIMERPATH : Filepath [DEFAULT: None]
        Filepath prefix for all primer outputs
    -------
    Converts MFEprimer dimer outputs into CSV tables
    """
    
    print("Reading in files............")
    # read in files - convert each interaction to a single line in array
    all_dimers = ReadDimerTXT(ALL_DIMERS)
    if END_DIMERS is None:
        end_dimers = []
    else:
        end_dimers = ReadDimerTXT(END_DIMERS)
    
    print("Remove duplicate dimers..................")
    # combine dimers into dataframe
    dimers_df = pd.concat([end_dimers, all_dimers])
    # get counts of each row
    counts = dimers_df.value_counts()
    # grab unique dimers between primers
    dimers = pd.DataFrame(list(counts.keys()), 
                          columns=['Primer1','Primer2','Pair1','Pair2','Score','DeltaG'])
    # grab unique sets of interacting primer pairs
    dimers_subset = pd.DataFrame(list(zip(dimers['Primer1'], dimers['Primer2'], dimers['Pair1'], dimers['Pair2'])),
                                 columns=['Primer1','Primer2', 'Pair1','Pair2'])
    counts_pairs = dimers_subset.value_counts()
    pair_dimers = pd.DataFrame(list(counts_pairs.keys()), columns=['Primer1','Primer2', 'Pair1','Pair2'])
    
    print("Extracting primer ID info.....")
    # get list of primer IDs, locusIDs, and primer pair IDs
    with open(ALL_DIMERS, 'r') as file:
        lines = file.readlines()
        start_indx = list(filter(lambda x: 'Primer ID' in lines[x], range(len(lines))))[0]
        end_indx = list(filter(lambda x: 'Dimer List' in lines[x], range(len(lines))))[0]
        primerIDs = [lines[x].split(' ')[0] for x in range(start_indx+3, end_indx-4)]
        #locusIDs = [primerIDs[x].split('.')[0]+'_'+primerIDs[x].split('.')[1] for x in range(len(primerIDs))]
        #pairIDs = [locusIDs[x]+'_'+primerIDs[x].split('_')[2] for x in range(len(primerIDs))]
    #del lines # clean up
    gc.collect() # clean up
    
    print("Calculating pairwise primer pair interactions........")
    pair_interactions = CalcPairwiseDimers(pair_dimers, primerIDs, pairs=True)
    
    # convert to binary format    
    pair_interactions_bin = pair_interactions.copy()
    pair_interactions_bin[pair_interactions_bin>=1] = 1
    
    print("Calculating total interactions per primer pair...........")
    # calculate sum of interactions for each primer pair
    pair_sums = pair_interactions.sum(axis=1)
    pair_sums_bin = pair_interactions_bin.sum(axis=1) # total pairs interacted with
    
    print("Saving primer pair output files!")
    ## Export all files
    # export total pairwise interactions per primer pair (wide)
    pairwide = OUTPATH + '_wide.csv'
    pair_interactions.to_csv(pairwide)
    
    # export binary pairwise interactions per primer pair (wide)
    pairwidebin = OUTPATH + '_binary_wide.csv'
    pair_interactions_bin.to_csv(pairwidebin)
    
    # export total interactions per primer pair (long)
    pairlong = OUTPATH + '_sum.csv'
    pair_sums.to_csv(pairlong)

    # export total # pairwise interactions per primer pair long (long)
    pairlongbin = OUTPATH + '_binary_sum.csv'
    pair_sums_bin.to_csv(pairlongbin)
    
    # calculate interactions per primer
    if OUTPRIMERPATH!="False":
        print("")
        print("")
        print("Calculating pairwise primer interactions..........")
        # primers - total # interactions
        dimers_sub = pd.DataFrame(list(zip(dimers['Primer1'],dimers['Primer2'],dimers['Pair1'],dimers['Pair2'])),
                                  columns=['Primer1','Primer2','Pair1','Pair2'])
        primer_interactions = CalcPairwiseDimers(dimers_sub, primerIDs, pairs=False)
    
        print("Converting primer interactions to binary format.......")
        primer_interactions_bin = primer_interactions.copy()
        primer_interactions_bin[primer_interactions_bin>=1] = 1
    
        # calculate # interactions per primer
        print("Calculating total interactions per primer............")
        primer_sums = primer_interactions.sum(axis=1) # total interactions       
        primer_sums_bin = primer_interactions.sum(axis=1) # total primers interacted with
        
        print("Saving primer output files!")
        # export total interactions per primer (long)
        primerlong = OUTPRIMERPATH + '_sum.csv'
        primer_sums.to_csv(primerlong)
        # export pairwise interactions per primer (wide)
        primerwide = OUTPRIMERPATH + '_wide.csv'
        primer_interactions.to_csv(primerwide)
        # export # primers each primer interacts with (long)
        primerlong = OUTPRIMERPATH + '_binary_sum.csv'
        primer_sums_bin.to_csv(primerlong)
        # export pairwise interactions per primer (wide)
        primerwide = OUTPRIMERPATH + '_binary_wide.csv'
        primer_interactions_bin.to_csv(primerwide)



def CalcPairwiseDimers(dimers, primerIDs, pairs=True):
    # duplicate dimers with pair1/primer1 and pair2/primer2 flipped (dimers will be counted by pair1 column, so 
    # this ensures that they'll be counted from both directions)
    flipped = pd.DataFrame(list(zip(dimers['Primer2'], dimers['Primer1'], dimers['Pair2'], dimers['Pair1'])),
                           columns=['Primer1','Primer2','Pair1','Pair2'])
    # combine raw and flipped dimers
    dimers_flipped = pd.concat([dimers, flipped])
    # add column for dimer count
    dimers_flipped['Count']=1
    # identify missing primers (i.e., primers with 0 dimers in files)
    primersRow = pd.unique(dimers_flipped['Primer1'])
    primersCol = pd.unique(dimers_flipped['Primer2'])
    missingRows = list(itertools.filterfalse(lambda x: x in primersRow, primerIDs))
    missingCols = list(itertools.filterfalse(lambda x: x in primersCol, primerIDs))
    # add rows for missing pairs (otherwise they won't be included in table)
    for primer in missingRows:
        pair = primer.replace('.REV', '').replace('.FW', '')
        dimers_flipped.loc[len(dimers_flipped)] = [primer, primer, pair, pair, 0]
    for primer in missingCols:
        pair = primer.replace('.REV', '').replace('.FW', '')
        dimers_flipped.loc[len(dimers_flipped)] = [primer, primer, pair, pair, 0]
    if pairs:
        # cross-tabulate pair1 vs pair2 counts
        pair_interactions = pd.crosstab(index=dimers_flipped['Pair1'], columns=dimers_flipped['Pair2'],
                                        dropna=False, values=dimers_flipped['Count'], aggfunc="sum")
    else:
        pair_interactions = pd.crosstab(index=dimers_flipped['Primer1'], columns=dimers_flipped['Primer2'],
                                        dropna=False, values=dimers_flipped['Count'], aggfunc="sum")
    # replace NAs with 0s
    pair_interactions = pair_interactions.fillna(0)
    # convert to integer
    pair_interactions = pair_interactions.astype(int)
    return pair_interactions



def ReadDimerTXT(infile):
    dimers = []
    with open(infile, 'r') as file:
        lines=file.readlines()
        dimer_indx = list(filter(lambda x: 'Dimer' in lines[x], range(len(lines))))[2:] # skip first two (headers)
        for i in dimer_indx:
            # grab all lines associated with this interaction & parse
            line1 = lines[i].strip().split(' ')
            line2 = lines[i+2].strip().split(' ')
            primer1 = line1[2]
            primer2 = line1[4]
            # we'll extract these later to reduce some memory reqs
            pair1 = primer1.replace('.FWD','').replace('.REV','')
            pair2 = primer2.replace('.FWD','').replace('.REV','')
            score = line2[1].replace(',','')
            delta_g = line2[5]
            #structure = lines[i+4:i+7] # no use for this yet
            dimers.append([primer1, primer2, pair1, pair2, score, delta_g])
    # convert to numpy array
    df = pd.DataFrame(dimers, columns=['Primer1','Primer2','Pair1','Pair2','Score','DeltaG'])
    return df
    gc.collect()#clean up



if __name__ == "__main__":
   main(sys.argv[1],
        sys.argv[2],
        sys.argv[3],
        str(sys.argv[4]))
    # set up pool specs for multiprocessing
   # with multiprocessing.Pool(4) as pool:
   #     print(Pooling.tabulateDimers)
