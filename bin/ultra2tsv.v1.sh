#!/usr/bin/env bash
#===============================================================================
#: TITLE:        json_to_tsv.sh
#: DESCRIPTION:  Convert JSON file(s) to TSV by extracting a fixed set of fields.
#: USAGE:        json_to_tsv.sh <json_file|glob>
#: EXAMPLES:     json_to_tsv.sh results/sample.json
#:               json_to_tsv.sh results/*.json
#: INPUTS:       One or more .json files. The script expects keys including
#:               SequenceName, Start, Length, Period, Score, Substitutions,
#:               Insertions, Deletions, Consensus, Sequence.
#: OUTPUTS:      For each input F.json, writes F.tsv in the same directory.
#: DEPENDENCIES: bash, perl, awk, coreutils
#: EXIT CODES:   0 success; 1 usage/help; 2 processing error
#: AUTHOR:       Brandon Jordan
#: LICENSE:      MIT
#: DATE:         2025-10-14
#===============================================================================

print_help() { grep -E '^#:' "$0" | sed 's/^#:[ ]\?//'; }

usage() {
    echo "Usage: $0 <json_file|glob>"
    echo "Try:   $0 --help"
    exit 1
}

# Check if no arguments or help requested
if [[ $# -eq 0 ]] || [[ $1 == "--help" ]] || [[ $1 == "-h" ]]; then
    print_help
    exit 1
fi

# processes the first argument (can be a glob pattern)
for file in "$@"; do
       	ref=`echo $file | perl -pe 's/(^*)\.json/$1.tsv/'` ;
       	cat $file | perl -pe 's/"//g; s/,//; s/\{//; s/\}//; s/://' | awk -v ORS="" '$1~/SequenceName/ {print "\n"} $1~/Start|Length|Period|Score|Substitutions|Insertions|Deletions|Consensus|Sequence/ {print $2 "\t"}' > $ref ; 
done
