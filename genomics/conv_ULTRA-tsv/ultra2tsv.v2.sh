#!/bin/bash

# Function to print usage instructions
usage() {
    echo "Usage: $0 <json_file>"
    echo "Converts a JSON file to a TSV file, extracting specific fields."
    exit 1
}

# Check if no arguments are provided or help is requested
if [[ $# -eq 0 ]] || [[ $1 == "--help" ]]; then
    usage
fi

for file in $1; do echo $file;
       	ref=`echo $file | perl -pe 's/(^*)\.json/$1.tsv/'` ;
       	cat $file | perl -pe 's/"//g; s/,//; s/\{//; s/\}//; s/://' | awk -v ORS="" '$1~/SequenceName/ {print "\n"} $1~/Start|Length|Period|Score|Substitutions|Insertions|Deletions|Consensus|Sequence/ {print $2 "\t"}' > $ref ; done
