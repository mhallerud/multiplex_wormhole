#!/bin/bash

# set up limits based on input params
Tm_limit=$1
dG_hairpins=$2
dG_end_limit=$3
dG_self_limit=$4
OUTDIR=$5
OUTNAME=$6


## filter acceptable primers based on secondary structures deltaG and annealing temps

# set up array to store primers that pass filtering
filtered_primers=()
filtered_primers+=("PrimerID,LocusID,PrimerPair,Direction,Sequence,StartBP,Length,AnnealTempC,PropBound,AmpliconSize") #CSV header
ids=() #array for locus IDs that pass filtering

# loop through each locus
for locus in "$(ls -1 $OUTDIR/1_InitialPrimers/*.out)"
do
	# progress tracking for # primers per locus
	passed=0

	# loop through each primer pair
	primerpairs=$(grep PRIMER_PAIR_NUM_RETURNED $locus | cut -d"=" -f 2)
	echo $primerpairs
	locusID=$(basename $locus | cut -d'.' -f 1)
	if [ "$primerpairs" -gt 0 ]; then
		echo "Testing " $locusID "........................................."
		for N in $(seq 0 $(($primerpairs-1)))
		do
			echo "     Primer #" $N
			
			# start counter to make sure all tests are passed
			declare -i tests=0
		
			# check FW and REV structures separately
			for DIR in LEFT RIGHT
			do
				if [ "$DIR" = "LEFT" ]; then declare DIRNAME="FW"; else declare DIRNAME="REV"; fi
			
				SelfAny=$(grep PRIMER_"$DIR"_"$N"_SELF_ANY_STUCT $locus)
				if [ "$SelfAny" != "" ]
				then
					echo "         " $DIRNAME " Self any secondary structure"
					# extract annealing temp of structure
					Tm="$(echo $SelfAny | cut -d"=" -f 2 | awk 'BEGIN{FS=";"}{print $1}' | head -n 1 | cut -d":" -f 2 | cut -d"&" -f 1)"
					Tm_int="$( printf "%.0f" $Tm )"
					# extract delta G of structures
					dG="$(echo $SelfAny | cut -d"=" -f 2 | awk 'BEGIN{FS=";"}{print $2}' | head -n 1 | cut -d" " -f 3)"
					dG_int="$( printf "%.0f" $dG)"
					# check if criteria are met
					#echo $dG_int $Tm_int
					if (($Tm_int > $Tm_limit)) && (($dG_int < $dG_self_limit))
					then
						echo "            FAILED"
					else
						tests+=1
					fi
				else
					tests+=1 #passed test bc this structure doesn't exist
				fi
			
				# check primer self end secondary structures
	       	      		SelfEnd=$(grep PRIMER_"$DIR"_"$N"_SELF_END_STUCT $locus)
	                	if [ "$SelfEnd" != "" ] #115
	                	then
	                       		echo "         " $DIRNAME "Self end secondary structure"
       	                 		# extract annealing temp of structure
                        		Tm="$(echo $SelfEnd | cut -d"=" -f 2 | awk 'BEGIN{FS=";"}{print $1}' | head -n 1 | cut -d":" -f 2 | cut -d"&" -f 1)"
					Tm_int="$( printf "%.0f" $Tm )"
					# extract delta G of structures
                        		dG="$(echo $SelfEnd | cut -d"=" -f 2 | awk 'BEGIN{FS=";"}{print $2}' | head -n 1 | cut -d" " -f 3)"
					dG_int="$( printf "%.0f" $dG )"					
					# check if criteria are met
                                        #echo $dG_int $Tm_int
                        		if (( $Tm_int > $Tm_limit)) && (($dG_int < $dG_end_limit))
                        		then
                                		echo "            FAILED"
                        		else
						tests+=1
                        		fi
				else
					tests+=1 #passed test bc this structure doesn't exist
				fi
                  
				# check primer hairpin structures
                		Hairpin=$(grep PRIMER_"$DIR"_"$N"_HAIRPIN_STUCT $locus)
                		if [ "$Hairpin" != "" ]
                		then
                        		echo "         " $DIRNAME "Hairpin structure"
                        		# extract annealing temp of structure
                        		Tm="$(echo $Hairpin | cut -d"=" -f 2 | awk 'BEGIN{FS=";"}{print $1}' | head -n 1 | cut -d":" -f 2 | cut -d"&" -f 1)"
                        		Tm_int="$( printf "%.0f" $Tm )"
					# extract delta G of structures
                        		dG="$(echo $Hairpin | cut -d"=" -f 2 | awk 'BEGIN{FS=";"}{print $2}' | head -n 1 | cut -d" " -f 3)"
					dG_int="$( printf "%.0f" $dG )"
					# check if criteria are met
                                        #echo $dG_int $Tm_int
                        		if (($Tm_int > $Tm_limit)) && (($dG_int < $dG_hairpins))
                        		then
                               	 		echo "            FAILED"
                       			else
						tests+=1
					fi
				else
					tests+=1 #passed test bc this structure doesn't exist
                        	fi
			done
			
			# check primer pair for secondary structures
			pairStruct=$(grep PRIMER_PAIR_"$N"_COMPL_ANY_STUCT $locus)
			if [ "$pairStruct" != "" ]
                	then
				echo "          Primer pair secondary structure"
				# extract annealing temp of structure
				Tm="$(echo $pairStruct | cut -d"=" -f 2 | awk 'BEGIN{FS=";"}{print $1}' | head -n 1 | cut -d":" -f 2 | cut -d"&" -f 1)"
				Tm_int="$( printf "%.0f" $Tm )"
				# extract delta G of structures
				dG="$(echo $pairStruct | cut -d"=" -f 2 | awk 'BEGIN{FS=";"}{print $2}' | head -n 1 | cut -d" " -f 3)"
				dG_int="$( printf "%.0f" $dG )"
				# check if criteria are met
                                #echo $dG_int $Tm_int
				if (($Tm_int > $Tm_limit)) && (($dG_int < $dG_self_limit))
				then    
					echo "            FAILED"
				else    
					tests+=1
				fi
			else
				tests+=1 #passed test bc this structure doesn't exist
			fi
		
			# check primer pair for end structures
	        	pairEndStruct=$(grep PRIMER_PAIR_"$N"_COMPL_END_STUCT $locus)
                	if [ "$pairEndStruct" != "" ]
                	then
                        	echo "          Primer pair end secondary structure"
                        	# extract annealing temp of structure
                        	Tm="$(echo $pairEndStruct | cut -d"=" -f 2 | awk 'BEGIN{FS=";"}{print $1}' | head -n 1 | cut -d":" -f 2 | cut -d"&" -f 1)"
				Tm_int="$( printf "%.0f" $Tm )"
				# extract delta G of structures
                        	dG="$(echo $pairEndStruct | cut -d"=" -f 2 | awk 'BEGIN{FS=";"}{print $2}' | head -n 1 | cut -d" " -f 3)"
                        	dG_int="$( printf "%.0f" $dG )"
				# check if criteria are met
                                #echo $dG_int $Tm_int
				if (($Tm_int > $Tm_limit)) && (($dG_int < $dG_end_limit))
                        	then
                        	        echo "            FAILED"
                        	else
                        	        tests+=1
                        	fi                     
                	else
                        	tests+=1 #passed test bc this structure doesn't exist
                	fi
		
			# if all tests were passed (n=8), then the primer sequences are moved into the filtered primers folder
			if [ $tests -eq 8 ]
			then
				# find elements that we want to save
				ID=$(basename $locus | cut -d"." -f 1)
				FWseq=$(grep PRIMER_LEFT_"$N"_SEQUENCE $locus | cut -d"=" -f 2 | tr '[:upper:]' '[:lower:]')
				FWpos=$(grep PRIMER_LEFT_"$N"= $locus | cut -d"=" -f 2)
				FWtm=$(grep PRIMER_LEFT_"$N"_TM $locus | cut -d"=" -f 2)
				FWbound=$(grep PRIMER_LEFT_"$N"_BOUND $locus | cut -d"=" -f 2)
				REVseq=$(grep PRIMER_RIGHT_"$N"_SEQUENCE $locus | cut -d"=" -f 2 | tr '[:upper:]' '[:lower:]')
				REVpos=$(grep PRIMER_RIGHT_"$N"= $locus | cut -d"=" -f 2)
				REVtm=$(grep PRIMER_RIGHT_"$N"_TM $locus | cut -d"=" -f 2)
				REVbound=$(grep PRIMER_RIGHT_"$N"_BOUND $locus | cut -d"=" -f 2)
				ampSize=$(grep PRIMER_PAIR_"$N"_PRODUCT_SIZE $locus | cut -d"=" -f 2)

				# add FW and REV primers separately to the array
				filtered_primers+=("$(echo $locusID"_"$N"_FW,"$locusID","$N",FW,"$FWseq","$FWpos","$FWtm","$FWbound","$ampSize)")
				filtered_primers+=("$(echo $locusID"_"$N"_REV,"$locusID","$N",REV,"$REVseq","$REVpos","$REVtm","$REVbound","$ampSize)")

				# progress tracking for the number of primers per locus that passed
				passed+=1
			fi

		done
	fi
	
	# add locus ID to list if any primers exist
	if (($passed > 0)); then
		ids+=($locusID)
	fi
done


## Export filtered primers as CSV
cd "$OUTDIR"
rm "$OUTNAME".csv
for i in "${!filtered_primers[@]}"
do
	printf "%s\n" "${filtered_primers[$i]}"
done > "$OUTNAME".csv


## Export locus IDs with primers as text file
rm "$OUTNAME"_LocusIDs.txt
for i in "${!ids[@]}"
do      
        printf "%s\n" "${ids[$i]}"
done > "$OUTNAME"_LocusIDs.txt
