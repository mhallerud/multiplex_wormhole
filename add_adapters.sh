#!/bin/bash

IN=$1
OUT=$2

# define a function to reverse complement DNA sequences
## via dacarlin/rc.bash https://gist.github.com/dacarlin/65965dfe130f0a2c4bbb
function ReverseComplement () {
        echo $1 | rev | tr atcg tagc
}

# get number of lines in primer CSV
lines="$(wc -l $IN | awk '{print $1}')"

# make empty array to store adapter-ligated primer sequences
adaptprime=()

# loop through each primer and add adapter
for i in $(seq 2 $lines) #line 1 is header
do
        # pull out row of CSV and primer direction
        DIR="$(awk -v i=$i 'NR==i' $IN | awk '{print $0}' | awk 'BEGIN{FS=","}{print $4}')"
	primerid="$(awk -v i=$i 'NR==i' $IN | awk '{print $0}' | awk 'BEGIN{FS=","}{print $1}')"
        
	# append adapter to primer sequence
        if [ "$DIR" = "FW" ]; then
		# extract primer from CSV
		primer="$(awk -v i=$i 'NR==i' $IN | awk '{print $0}' | awk 'BEGIN{FS=","}{print $5}')"
                # add adapter to primer
		primerAdapt=$(echo "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG""$primer")
                # convert to all lowercase
		primerAdapt="$(echo $primerAdapt | tr '[:upper:]' '[:lower:]')"
		#echo $primerAdapt # progress checking / error checking
		# add primer ID to list (this will make the output in fasta format)
		adaptprime+=($(echo ">"$primerid))
		# add adapter-ligated primer sequence to list
		adaptprime+=($primerAdapt)
        fi
	
	# do the same for reverse primers, except that the primer sequence gets reverse-complemented prior to adding the adapter        
        if [ "$DIR" = "REV" ]; then
                primer="$(echo $(awk -v i=$i 'NR==i' $IN | awk '{print $0}' | awk 'BEGIN{FS=","}{print $5}'))"
		primer_rc="$(ReverseComplement $primer)"
		primerAdapt_rc=$(echo "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG""$primer_rc")
                primerAdapt_rc="$(echo $primerAdapt_rc | tr '[:upper:]' '[:lower:]')" #change to all lowercase
		#echo $primerAdapt_rc
		adaptprime+=($(echo ">"$primerid))
		adaptprime+=($primerAdapt_rc)
        fi
done
 
# export adapter-primer sequences to CSV
rm $OUT
for i in "${!adaptprime[@]}"
do
        printf "%s\n" "${adaptprime[$i]}"
done > $OUT

# save primers and IDs separately
grep ">" $OUT > "$OUT"_primerIDs.temp
grep -v ">" $OUT > "$OUT"_primerSeqs.temp
