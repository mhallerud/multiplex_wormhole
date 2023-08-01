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
./scripts/primer3_batch_design.sh "$TEMPLATES" "$OUTDIR" > primer3_batch.log
numMissing="$(grep 'No primers found' primer3_batch.log | wc -l | awk '{print $1}')"
numWithPrimers=$(($numIn-$numMissing))
echo "      # loci primers could be designed for: " $numWithPrimers


# filter primers
## dimer melting temp=Tm < 45 Celsius (min primer Tm - 10)
## deltaG <= -3 kcal/mol for any hairpins
## deltaG <= -5 kcal/mol for any dimers at the 3' ends (pair or self)
## deltaG <= -10 kcal/mole in the middle of primers where steric hindrance makes structures difficult to form - also for pair heterodimers in middle
echo " "
echo "Filtering primers based on delta G and melting temp (Tm)........" 
./scripts/filter_primers_Tm_dG.sh 45 -3000 -5000 -10000 "$OUTDIR" "$NAME"_filteredPrimers > filter_primers.log
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
./scripts/convert_CSV_to_fasta.sh "$OUTDIR"/"$NAME"_specificityCheck_passed.csv 1 5
echo " " 
echo "NEXT STEP: Input to PrimerSuite PrimerDimer (http://www.primer-dimer.com/): " "$OUTDIR"/"$NAME"_specificityCheck_passed.csv


# check for primer dimers using SADDLE & PrimerSuite PrimerDimer function (http://www.primer-dimer.com/)
# create symbolic link for matlab & run SADDLE
#ln -s /Applications/MATLAB_R2022b.app/bin/matlab matlab
# input $NAME_primerAdapters.fa to PrimerSuite
