#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 08:51:53 2023

@author: maggiehallerud
"""

IN='/Users/maggiehallerud/Marten_Primer_Design/TEST/PrimerDimerReport_25jul2023.txt'

# grab out data for each interaction
primer1_id = []
primer2_id = []
deltaG = []
primer1_seq = []
bonds = []
with open(IN, "r") as file:
    lines = file.readlines()
    for line in lines:
        if 'kcal/mol' in line:
            linesplit=line.split(' vs ')
            primer1=linesplit[0].strip()
            primer2=linesplit[1].split(': ')[0].strip()
            deltaG_=float(line.split(': ')[1].split(' kcal/mol')[0])
            primer1_id.append(primer1)
            primer2_id.append(primer2)
            deltaG.append(deltaG_)
        if "5'>" in line:
            primer1_seq.append(line.strip())
        if "|" in line:
            bonds.append(line)

# convert bonds line to bp position for each bond, then
# extract nucleotides bonded based on position
bond_pos=[]
bond_nucs=[]
GC = 
for row in range(len(bonds)):
    # extract bp position for each bond
    bond = bonds[row]
    pos_indx = []
    bondsplit = bond.split(' ')[11:]#skip past '      5'> ' prefix
    bond_indx = list(filter(lambda x: '|' in bondsplit[x], range(len(bondsplit))))
    i=1
    for b in bond_indx:
        n_bonds=len(bondsplit[b])
        b_inds = list(range(b+i,b+i+n_bonds))
        for ind in b_inds:#append indices individually
            pos_indx.append(ind)# this is a list of all bp positions bonded
        i+=n_bonds#one space is lost between every set of bonds
    bond_pos.append(pos_indx)
    # extract nucleotide and sequence info for each
    seq = primer1_seq[row].split('>')[1].strip()
    bond_nucs.append([seq[i] for i in pos_indx])
    bond_pos.append(pos_indx)

# TODO: Add melting temperature for primer dimers: Tm = 81.5 + 0.41(%GC) - 675/N - % mismatch

## Calculate melting temperature for primer dimers: Tm = 81.5 + 0.41(%GC) - 675/N - % mismatch


# grab all primer interactions
for (line in 1:length(dimers)){
  
}#line
