#!/bin/bash

# This script will discard any primer pairs that are going to amplify multiple SNPs

IN=$1
OUT=$2
TEMPLATES=$3

# define a function to reverse complement DNA sequences
## via dacarlin/rc.bash https://gist.github.com/dacarlin/65965dfe130f0a2c4bbb
function ReverseComplement () {
	echo $1 | rev | tr atcg tagc
}


# get number of lines in primer CSV
lines="$(wc -l $IN | awk '{print $1}')"

# check that each primer only matches to 1 line in templates
failedIDs=() #array that will hold failed primer pair IDs
failedSeqs=() #array that will hold failed primer sequences

for i in $(seq 2 $lines)
do
	# find lines including that primer (adapters trimmed)
	primer="$(awk -v i=$i 'NR==i' $IN | awk '{print $0}' | awk 'BEGIN{FS=","}{print $5}' | 
sed 's/tcgtcggcagcgtcagatgtgtataagagacag//' | sed 's/gtctcgtgggctcggagatgtgtataagagacag//')"
	# get the number of matching lines
	n_matches="$(grep $primer $TEMPLATES | wc -l | awk '{print $1}')"
	
	# if the number of matches is 0, try the reverse complement
	if (($n_matches==0))
	then
		primer_rc="$(ReverseComplement $primer)"
		n_matches="$(grep $primer_rc $TEMPLATES | wc -l | awk '{print $1}')"
	fi

	# if there's more than 1 match, the primer can't be used because it's not specific to a given SNP
	if (($n_matches>1)); then
		primerPairID="$(awk -v i=$i 'NR==i' $IN | awk '{print $0}' | awk 'BEGIN{FS=","}{gsub("_FW", ""); gsub("_REV",""); print $1}')"
		failedIDs+=($primerPairID)
		failedSeqs+=($primer)
	fi
done

# export failed primer sequences to TXT
rm "$OUT"_failedIDs.txt
for i in "${!failedIDs[@]}"
do
        printf "%s\n" "${failedIDs[$i]}"
done > "$OUT"_failedIDs.txt

# export failed primer pair IDs to TXT
rm "$OUT"_failedSeqs.txt
for i in "${!failedSeqs[@]}"
do
        printf "%s\n" "${failedSeqs[$i]}"
done > "$OUT"_failedSeqs.txt

# save CSV for primers that passed
rm "$OUT"_passed.csv
grep -v -f "$OUT"_failedIDs.txt $IN > "$OUT"_passed.csv


# save FASTA for primers that passed
lines="$(wc -l "$OUT"_passed.csv | awk '{print $1}')"
passed=()
for i in $(seq 2 $lines)
do
	primer="$(awk -v i=$i 'NR==i' "$OUT"_passed.csv | awk 'BEGIN{FS=","}{print $5}')"
	id="$(awk -v i=$i 'NR==i' "$OUT"_passed.csv | awk 'BEGIN{FS=","}{print $1}')"
	passed+=("$(echo \>$id)")
	passed+=($primer)
done

rm "$OUT"_passed.fa
for i in "${!passed[@]}"
do
        printf "%s\n" "${passed[$i]}"
done > "$OUT"_passed.fa