#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 11:57:09 2023

@author: maggiehallerud
"""

# load dependencies
import os
import sys

# import user-defined parameters from command line
ALL_DIMERS = sys.argv[1]
END_DIMERS = sys.argv[2]
OUTNAME = sys.argv[3]
ALL_DIMERS='/Users/maggiehallerud/Marten_Primer_Design/examples/Example_MFEprimerDimers.txt'
END_DIMERS='/Users/maggiehallerud/Marten_Primer_Design/examples/Example_MFEprimerDimers_ends.txt'
OUTNAME='Example'



def main():
    # read in files - convert each interaction to a single line in array
    all_dimers = readDimerFile(ALL_DIMERS)
    if END_DIMERS == '':
        end_dimers = []
    else:
        end_dimers = readDimerFile(END_DIMERS)
    
    # subset all dimers for all non-end dimers
    middle_indx = list(filter(lambda x: all_dimers[x] not in end_dimers, range(len(all_dimers))))
    middle_dimers = [all_dimers[x] for x in middle_indx]
    
    # combine middle and end dimers into same array
    dimers = end_dimers + middle_dimers
    
    # get list of primer IDs, locusIDs, and primer pair IDs
    with open(ALL_DIMERS, 'r') as file:
        lines = file.readlines()
        start_indx = list(filter(lambda x: 'Primer ID' in lines[x], range(len(lines))))[0]
        end_indx = list(filter(lambda x: 'Dimer List' in lines[x], range(len(lines))))[0]
        primerIDs = [lines[x].split(' ')[0] for x in range(start_indx+3, end_indx-4)]
        locusIDs = [primerIDs[x].split('_')[0]+'_'+primerIDs[x].split('_')[1] for x in range(len(primerIDs))]
        pairIDs = [locusIDs[x]+'_'+primerIDs[x].split('_')[2] for x in range(len(primerIDs))]

    # convert line by line dimers into table
    # total # interactions
    primer_interactions, primer_interactions_bin = tabulateDimers(dimers, primerIDs, locusIDs, pairs=False)
            
    # aggregate the number of interactions per primer pair
    primer_pairs = list(set(pairIDs))
    pair_loci = [primer_pairs[x].split('_')[0]+'_'+primer_pairs[x].split('_')[1] for x in range(len(primer_pairs))]
    pair_interactions, pair_interactions_bin = tabulateDimers(dimers, primer_pairs, pair_loci, pairs=True)        

    # calculate # interactions per primer (long format)
    #skipped- add later if it seems useful
    
    # calculate # interactions per primer pair (long format)
    pair_sums = totalDimers(pair_interactions, primer_pairs)
        
    # calculate # primer pairs each pair interacts with (binary long format)
    pair_bin_sums = totalDimers(pair_interactions_bin, primer_pairs)
    
    # Export all files
    OUTDIR = os.path.dirname(ALL_DIMERS)

    # export total pairwise interactions per primer pair (wide)
    pairwide = os.path.join(OUTDIR, OUTNAME+'_PrimerPairInteractions_wide.csv')
    exportToCSV(pair_interactions, pairwide)
    
    # export binary pairwise interactions per primer pair (wide)
    pairwidebin = os.path.join(OUTDIR, OUTNAME+'_PrimerPairInteractions_wide_binary.csv')
    exportToCSV(pair_interactions_bin, pairwidebin)
    
    # export total interactions per primer pair (long)
    pairlong = os.path.join(OUTDIR, OUTNAME+'_PrimerPairInteractions_sum.csv')
    exportToCSV(pair_sums, pairlong)

    # export total # pairwise interactions per primer pair long (long)
    pairlongbin = os.path.join(OUTDIR, OUTNAME+'_PrimerPairInteractions_sum_binary.csv')
    exportToCSV(pair_bin_sums, pairlongbin)
    
    ## export total interactions per primer (long)
    #primerlong = os.path.join(OUTDIR, OUTNAME+'_PrimerInteractions_sum.csv')
    #exportToCSV(primer_interactions_sum, primerlong) #primer_interactions_sum would need to be calculated
    # skipped
    
    # export pairwise interactions per primer (wide)
    #primerwide = os.path.join(OUTDIR, OUTNAME+'_PrimerInteractions_wide.csv')
    #exportToCSV(primer_interactions, primerwide)
    #skipped
    


def readDimerFile(infile):
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
            pair1 = primer1.replace('_FW','').replace('_REV','')
            pair2 = primer2.replace('_FW','').replace('_REV','')
            score = line2[1].replace(',','')
            delta_g = line2[5]
            #structure = lines[i+4:i+7] # no use for this yet
            dimers.append([primer1, primer2, pair1, pair2, score, delta_g])
        return dimers



def tabulateDimers(dimers, IDlist, locusIDs, pairs=False):
    # set up arrays to hold counts and binary interactions
    pairwise_interactions = [["NA" for _ in range(len(IDlist)+1)] for _ in range(len(IDlist)+1)]
    pairwise_interactions[0] = ['Primer2_ID'] + IDlist
    pairwise_interactions_bin = pairwise_interactions
    # loop through each primer or pair
    for i in range(len(IDlist)):
        # the first index will be the 'rows'
        rowID = IDlist[i]
        pairwise_interactions[i+1][0] = pairwise_interactions_bin[i+1][0] = rowID
        # compare with every other primer/pairs
        for j in range(len(IDlist)):
            # these will be the 'columns'
            colID = IDlist[j]
            # if these primers are for the same locus, then set to 0 because
            # we already filtered for homodimers and pair dimers, and dimers 
            # between pairs for the same locus don't matter because there will 
            # only ever be one primer pair per locus in a set
            if locusIDs[i]==locusIDs[j]:
                pairwise_interactions[i+1][j+1] = 0 # offset indices by 1 to allow space for row and column names
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
                if Ndimers>0:
                    Ndimers_bin = 1
                else:
                    Ndimers_bin = 0
                pairwise_interactions[i+1][j+1] = Ndimers
                pairwise_interactions_bin[i+1][j+1] = Ndimers_bin
    return pairwise_interactions, pairwise_interactions_bin



def totalDimers(pairwise_interactions, ids):
    sums = []
    for id in ids:
        # grab the row for this id
        indx = list(filter(lambda x: pairwise_interactions[x][0]==id, range(len(pairwise_interactions))))
        dimers = [pairwise_interactions[x] for x in indx]
        # sum the dimer count across all pairwise interactions
        Ndimers = sum(dimers[1:]) # skip first record because it's the ID field
        # add to output array
        sums.append([id, Ndimers])
    return sums



def exportToCSV(inArray, outCSV):
    with open(outCSV, 'w') as file:
        for line in inArray:
            str_line = str(line)[1:-1].replace("'", "")
            file.write(str_line+"\n")



if __name__ == "__main__":
    main()
