#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 17:43:26 2023

@author: maggiehallerud
"""
import multiprocessing
import sys



def main(i, primer_ids, locus_ids, dimers, pairs):
    # set up emty table
    pairwise_interactions = []
    pairwise_interactions.append(primer_ids)    
    return i, primer_ids, locus_ids, dimers, pairs, pairwise_interactions



def tabulateByRow(i, primer_ids, locus_ids, dimers, pairs):
    # set up empty arrays to hold results
    interactions_row = []
    # add rowname as first value in array
    rowID = primer_ids[i]
    interactions_row.append(rowID)
    # loop through every other primer pair to find dimers- these will be the 'columns'
    for j in range(len(primer_ids)):
        # grab primer ID
        colID = primer_ids[j]
        # if these primers are for the same locus, then set to 0 because
        # we already filtered for homodimers and pair dimers, and dimers 
        # between pairs for the same locus don't matter because there will 
        # only ever be one primer pair per locus in a set
        if locus_ids[i]==locus_ids[j]:
            interactions_row.append(0)
        else:
            # for pairs, look at fields 3-4 in dimer array
            if pairs:
                # get all dimers for these primers (including both pairwise comparisons)
                sub1_dimer_indx = list(filter(lambda x: dimers[x][2]==rowID and dimers[x][3]==colID, range(len(dimers))))
                sub2_dimer_indx = list(filter(lambda x: dimers[x][2]==rowID and dimers[x][3]==colID, range(len(dimers))))
                sub_dimer_indx = sub1_dimer_indx + sub2_dimer_indx
            # for primers, look at fields 1-2 in dimer array
            else:
                sub1_dimer_indx = list(filter(lambda x: dimers[x][0]==rowID and dimers[x][1]==colID, range(len(dimers))))
                sub2_dimer_indx = list(filter(lambda x: dimers[x][0]==rowID and dimers[x][1]==colID, range(len(dimers))))
                sub_dimer_indx = sub1_dimer_indx + sub2_dimer_indx
            # grab the # of primer interactions for this comparison
            Ndimers = len(sub_dimer_indx)
            # add to row
            interactions_row.append(Ndimers)
    # return populated row
    return(interactions_row)



if __name__=="__main__":
    # read in arguments from command line
    i = sys.argv[0] 
    primer_ids = sys.argv[1]
    locus_ids = sys.argv[2]
    dimers = sys.argv[3]
    pairs = sys.argv[4]
    # define arguments
    i, primer_ids, locus_ids, dimers, pairs, pairwise_interactions = main(i, primer_ids, locus_ids, dimers, pairs)
    # run tabulateByRow for each primer using multiprocessing pool
    with multiprocessing.Pool(4) as pool:
        results = [pool.apply(tabulateByRow, args=(x, primer_ids, locus_ids, dimers, pairs)) for x in range(len(primer_ids))]
        pairwise_interactions.append(results)
        print(pairwise_interactions)