IN=$1
ID_FIELD=$2
SEQ_FIELD=$3

OUTNAME="$(basename $IN | cut -d'.' -f 1)"
OUTDIR="$(dirname $IN)"
OUT="$(echo $OUTDIR/$OUTNAME)"

# get number of lines in input CSV
lines="$(wc -l $IN | awk '{print $1}')"

# make empty arrays to store ids and sequences
outlines=()

# loop through each primer and add adapter
for i in $(seq 2 $lines) #line 1 is header
do
	# grab sequences and IDs from CSV
	id="$(awk -v i=$i '(NR==i)' $IN | awk -v id=$ID_FIELD 'BEGIN{FS=","}{print $id}')"
	seq="$(awk -v i=$i '(NR==i)' $IN | awk -v seq=$SEQ_FIELD 'BEGIN{FS=","}{print $seq}')"
	# reformat id
	id2="$(echo '>'$id)"
	# append to lists
	outlines+=($id2)
	outlines+=($seq)
done

# write IDs and sequences to fasta file
for i in "${!outlines[@]}"
do
        printf "%s\n" "${outlines[$i]}"
done > "$OUT".fa


# save primers and IDs separately
grep ">" "$OUT".fa > "$OUT"_primerIDs.temp
grep -v ">" "$OUT".fa > "$OUT"_primerSeqs.temp
