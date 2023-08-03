#!/bin/bash

TEMPLATES=$1 #template sequences that primers are designed from
OUTDIR=$2 #directory where all output files will be stored
NAME=$3 #prefix for naming outputs

numIn="$(wc -l $TEMPLATES | awk '{print $1-1}')"
echo "# loci in input file: " $numIn

# design initial primers using Primer3
# NOTE: adapters are added to primers and considered in secondary structure calculations via the LEFT_PRIMER_OVERHANG & RIGHT_PRIMER_OVERHANG settings
# NOTE: the primer3 setting files are setup so that no filtering occurs based on secondary structures at this step
echo " "
echo "Designing initial primers...."
echo "      (see primer3_batch.log for details)"
./scripts/primer3_batch_design.sh "$TEMPLATES" "$OUTDIR" > "$OUTDIR"/primer3_batch.log
numMissing="$(grep 'No primers found' "$OUTDIR"/primer3_batch.log | wc -l | awk '{print $1}')"
numWithPrimers=$(($numIn-$numMissing))
echo "      # loci primers could be designed for: " $numWithPrimers


# filter primers
## dimer melting temp=Tm < 45 Celsius (min primer Tm - 10)
## deltaG <= -2 kcal/mol for any hairpins
## deltaG <= -5 kcal/mol for any dimers at the 3' ends (pair or self)
## deltaG <= -10 kcal/mole in the middle of primers where steric hindrance makes structures difficult to form - also for pair heterodimers in middle
echo " "
echo "Filtering primers based on delta G and melting temp (Tm)........" 
# in order, the arguments for this script are:
# Tm_limit : Max melting temp allowed for secondary structures
# #dG_hairpins : Min delta G allowed for hairpins
# dG_end_limit : Min delta G allowed for secondary structures at the end of primers (self-dimers or heterodimers)
# dG_self_limit : Min delta G allowed for secondary structures occurring in the middle of primers
# OUTDIR : Directory to save outputs to
# OUTNAME : Prefix for naming outputs
./scripts/filter_primers_Tm_dG.sh 45 -2000 -5000 -10000 "$OUTDIR" "$NAME"_filteredPrimers > "$OUTDIR"/filter_primers.log
numfiltered="$(wc -l "$OUTDIR"/"$NAME"_filteredPrimers_LocusIDs.txt | awk '{print $1}')"
echo "      # loci with primers after filtering:" $numfiltered


# check primer specificity, remove any non-specific primer pairs
echo " "
echo "Checking specificity of remaining primers......." 
./scripts/check_specificity.sh "$OUTDIR"/"$NAME"_filteredPrimers.csv "$OUTDIR"/"$NAME"_specificityCheck "$TEMPLATES"
echo "      Number of primer pairs that failed specificity check: " "$(sort "$OUTDIR"/"$NAME"_specificityCheck_failedIDs.txt | wc -l | awk '{print $1}')"
echo "      Number of unique primer sequences that failed specificity check: " "$(sort "$OUTDIR"/"$NAME"_specificityCheck_failedSeqs.txt | uniq | wc -l | awk '{print $1}')"
echo "      Loci represented in non-specific primers: " 
echo "$(cat "$OUTDIR"/"$NAME"_specificityCheck_failedIDs.txt| tr "_" " " | awk '{print $1"_"$2}' | sort | uniq)"


# add adapters - this is already done in primer3 with the PRIMER_*_OVERHANG settings!
#cd $OUTDIR
#./scripts/add_adapters.sh $NAME_filteredPrimers.csv $NAME_primerAdapters.fa


# convert checkSpecificity CSV to fasta format
#./scripts/convert_CSV_to_fasta.sh "$OUTDIR"/"$NAME"_specificityCheck_passed.csv 1 5
#echo " " 
#echo "NEXT STEP: Input to PrimerSuite PrimerDimer (http://www.primer-dimer.com/): " "$OUTDIR"/"$NAME"_specificityCheck_passed.csv


# originally, primers were checked via the PrimerSuite PrimerDimer function (http://www.primer-dimer.com/)
# PrimerSuite PrimerDimerReport files can be converted to the necessary table/sum files using scripts/translate_primerSuite_report.R

# use MFEprimer to calculate primer dimers
echo " " 
echo "Calculating dimers with MFEprimer dimer........."
# -i = input FASTA of primer sequences 
# -o = output file
# -d = maximum deltaG threshold to consider dimers (kcal/mol)
# -s = minimum score threshold to consider dimers(scores are calculated with +1 for each match and -1 for each mismatch (not including Ns)
# -m = max allowed mismatches per dimer
# -p = only output dimers with 3' end bind
# --diva = concentration of divalent cations (mM)
# --mono = concentration of monovalent cations (mM)
# --dntp = concentration of dNTPs (mM)
# --oligo = concentration of annealing oligos (nM) 
../mfeprimer-3.2.7-darwin-10.6-amd64 dimer -i "$OUTDIR"/"$NAME"_specificityCheck_passed.fa -o "$OUTDIR"/"$NAME"_MFEprimerDimers.txt -d -6 -s 3 -m 40 --diva 3.8 --mono 50 --dntp 0.25 --oligo 50
../mfeprimer-3.2.7-darwin-10.6-amd64 dimer -i "$OUTDIR"/"$NAME"_specificityCheck_passed.fa -o "$OUTDIR"/"$NAME"_MFEprimerDimers_ends.txt -d -2.5 -s 3 -m 40 --diva 3.8 --mono 50 --dntp 0.25 --oligo 50 -p

# convert MFEprimer dimer report to tables
echo " " 
echo "Converting dimer report to tables........."
#script arguments order: 1- path to all dimers file, 2- path to end dimers file, 3- output name
python ./scripts/tabulate_MFEprimer_dimers.py "$OUTDIR"/"$NAME"_MFEprimerDimers.txt "$OUTDIR"/"$NAME"_MFEprimerDimers_ends.txt "$NAME"

# run optimization
echo " "
echo "Running optimization process......."
OPTOUT="$OUTDIR"/2_OptimizedSets
mkdir $OPTOUT

# script arguments in order:
#1- PRIMER_FASTA : path to FASTA of adapter-ligated primer sequences
#2- DIMER_SUMS : path to dimer load totals for each primer (long format)
#3- DIMER_TABLE : path to pairwise dimer table
#4- N_LOCI : Number of SNPs desired in final panel
#5- ITERATIONS : Number of iterations
#6- RUN : Run ID

echo "......Run 1"
python ./scripts/optimize_primers.py "$OUTDIR"/"$NAME"_specificityCheck_passed.fa "$OUTDIR"/"$NAME"_PrimerPairInteractions_sum_binary.csv "$OUTDIR"/"$NAME"_PrimerPairInteractions_wide_binary.csv 50 5000 1 > "$OPTOUT"/"$NAME"_OptimizedSet1.txt
echo " "

echo "......Run 2"
python ./scripts/optimize_primers.py "$OUTDIR"/"$NAME"_specificityCheck_passed.fa "$OUTDIR"/"$NAME"_PrimerPairInteractions_sum_binary.csv "$OUTDIR"/"$NAME"_PrimerPairInteractions_wide_binary.csv 100 5000 2 > "$OPTOUT"/"$NAME"_OptimizedSet2.txt
echo " "
 
echo "......Run 3"
python ./scripts/optimize_primers.py "$OUTDIR"/"$NAME"_specificityCheck_passed.fa "$OUTDIR"/"$NAME"_PrimerPairInteractions_sum_binary.csv "$OUTDIR"/"$NAME"_PrimerPairInteractions_wide_binary.csv 100 5000 3 > "$OPTOUT"/"$NAME"_OptimizedSet3.txt
echo " "

echo "......Run 4"
python ./scripts/optimize_primers.py "$OUTDIR"/"$NAME"_specificityCheck_passed.fa "$OUTDIR"/"$NAME"_PrimerPairInteractions_sum_binary.csv "$OUTDIR"/"$NAME"_PrimerPairInteractions_wide_binary.csv 100 5000 4 > "$OPTOUT"/"$NAME"_OptimizedSet4.txt
echo " "

echo "......Run 5"
python ./scripts/optimize_primers.py "$OUTDIR"/"$NAME"_specificityCheck_passed.fa "$OUTDIR"/"$NAME"_PrimerPairInteractions_sum_binary.csv "$OUTDIR"/"$NAME"_PrimerPairInteractions_wide_binary.csv 100 5000 5 > "$OPTOUT"/"$NAME"_OptimizedSet5.txt
echo " "

echo "......Run 6"
python ./scripts/optimize_primers.py "$OUTDIR"/"$NAME"_specificityCheck_passed.fa "$OUTDIR"/"$NAME"_PrimerPairInteractions_sum_binary.csv "$OUTDIR"/"$NAME"_PrimerPairInteractions_wide_binary.csv 100 5000 6 > "$OPTOUT"/"$NAME"_OptimizedSet6.txt
echo " "

echo "......Run 7"
python ./scripts/optimize_primers.py "$OUTDIR"/"$NAME"_specificityCheck_passed.fa "$OUTDIR"/"$NAME"_PrimerPairInteractions_sum_binary.csv "$OUTDIR"/"$NAME"_PrimerPairInteractions_wide_binary.csv 100 5000 7 > "$OPTOUT"/"$NAME"_OptimizedSet7.txt
echo " "

echo "......Run 8"
python ./scripts/optimize_primers.py "$OUTDIR"/"$NAME"_specificityCheck_passed.fa "$OUTDIR"/"$NAME"_PrimerPairInteractions_sum_binary.csv "$OUTDIR"/"$NAME"_PrimerPairInteractions_wide_binary.csv 100 5000 8 > "$OPTOUT"/"$NAME"_OptimizedSet8.txt
echo " "

echo "......Run 9"
python ./scripts/optimize_primers.py "$OUTDIR"/"$NAME"_specificityCheck_passed.fa "$OUTDIR"/"$NAME"_PrimerPairInteractions_sum_binary.csv "$OUTDIR"/"$NAME"_PrimerPairInteractions_wide_binary.csv 100 5000 9 > "$OPTOUT"/"$NAME"_OptimizedSet9.txt
echo " "

echo "......Run 10"
python ./scripts/optimize_primers.py "$OUTDIR"/"$NAME"_specificityCheck_passed.fa "$OUTDIR"/"$NAME"_PrimerPairInteractions_sum_binary.csv "$OUTDIR"/"$NAME"_PrimerPairInteractions_wide_binary.csv 100 5000 10 > "$OPTOUT"/"$NAME"_OptimizedSet10.txt
echo " "
