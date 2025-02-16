#!/bin/bash
# Author: Brandon D. Jordan
# 2023
# Obj: To convert JSON output from ULTRA tandem repeat finder to a tab-delimited format. This program can accept a single or batch of JSON files.
for file in $1; do echo $file;
       	ref=`echo $file | perl -pe 's/(^*)\.json/$1.tsv/'` ;
       	cat $file | perl -pe 's/"//g; s/,//; s/\{//; s/\}//; s/://' |       awk -v ORS="" '$1~/SequenceName/ {print "\n"} $1~/Start|Length|Period|Score|Substitutions|Insertions|Deletions|Consensus|Sequence/ {print $2 "\t"}' > $ref
done
