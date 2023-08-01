#!/bin/bash

currdir="$(pwd)"
settings=$4
outdir=$5

# set up settings file
rm primer3_settings.txt

sed 's/SEQUENCE_ID=/SEQUENCE_ID='"$1"'/' $settings | \
sed 's/SEQUENCE_TEMPLATE=/SEQUENCE_TEMPLATE='"$2"'/' | \
sed 's/SEQUENCE_TARGET=/SEQUENCE_TARGET='"$3"'/' > primer3_settings.txt

# design primers with primer3
#primer3_core [ --format_output ] [--default_version=1|--default_version=2] [ --io_version=4 ] [ --p3_settings_file=<file_path> ] [ --echo_settings_file ] [ --strict_tags ] [ --output=<file_path> ] [ --error=<file_path> ] [ input_file ] 
cd $outdir
/Users/maggiehallerud/primer3/src/primer3_core "$currdir"/primer3_settings.txt --error=$1.err --output=$1.out
