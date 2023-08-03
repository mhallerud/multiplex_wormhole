#!/bin/bash

## set up output directory structure
[[ -z "$1" ]] && { echo "Please provide an output directory:"; echo "primer3_batch_design.sh <OUTDIR>" ; exit 1; }

TEMPLATES=$1
OUTDIR=$2
#outls="$(ls -1 $OUTDIR | wc -l | awk '{print $1}')"
#if (($outls == 0))
#then
#	mkdir $OUTDIR
#	OUT1=$OUTDIR/1_InitialPrimers
#	mkdir $OUT1
#fi

OUT1=$OUTDIR/1_InitialPrimers
mkdir $OUT1
rm $OUT1/*
 
## load CSV file with SEQUENCE_ID, SEQUENCE_TEMPLATE, SEQUENCE_TARGET in columns
# use create_in_templates.R script lines 1-90 to create this CSV
# read in CSV
echo "Reading in sequences"
arr_csv=()
declare -i index=0
while IFS= read -r line
do
	arr_csv+=("$line")
	index+=1
done < $TEMPLATES


## design FW and REV primers for each sequence
echo "Designing initial primers"
index="$((index - 1))" #not sure why, but there's an extra line in arr_csv. This fixes it.
for row in $(seq 2 $index) #start at 2 to skip header
do
	# pull in inputs
	id=$(echo ${arr_csv[$row]} | awk 'BEGIN{FS=","}{gsub(/"/,""); print $1}')
	seq=$(echo ${arr_csv[$row]} | awk 'BEGIN{FS=","}{gsub(/"/,""); print $2}')
	SNP=$(echo ${arr_csv[$row]} | awk 'BEGIN{FS=","}{gsub(/"/,""); print $3","$4}')

	# run primer3 based on strict settings
	echo "-Designing primers for " $id
	./scripts/primer3.sh "$id" "$seq" "$SNP" ./primer3_settings/primer3_Base_NoSecondaryFilters.txt $OUT1

	# rerun with broad settings if no primers were found
	primers="$(grep PRIMER_PAIR_NUM_RETURNED=0 $OUT1/$id.out)"
	if [ "$primers" != "" ]
	then
		echo ".....Retrying with broader parameters"
		./scripts/primer3.sh $id $seq $SNP ./primer3_settings/primer3_Broad_NoSecondaryFilters.txt $OUT1
	fi

	# raise message if no primers were found for broad settings either
	primers="$(grep PRIMER_PAIR_NUM_RETURNED=0 $OUT1/$id.out)"
	if [ "$primers" != "" ]
	then
		echo ".....No primers found"
	fi

	# raise message if there were errors in primer design
	errors="$(wc -l $OUT1/$id.err | awk '{print $1}')"
	outErrors="$(grep ERROR $OUT1/$id.out)"
	if [ $errors  -gt 0 ] || [ "$outErrors" != "" ]
	then
		echo ".....ERROR in primer design!"
	fi
done
