#!/bin/bash

for file in $1; do echo $file;
       	ref=`echo $file | perl -pe 's/(^*)\.json/$1.tsv/'` ;
       	cat $file | perl -pe 's/"//g; s/,//; s/\{//; s/\}//; s/://' | awk -v ORS="" '$1~/SequenceName/ {print "\n"} $1~/Start|Length|Period|Score|Substitutions|Insertions|Deletions|Consensus|Sequence/ {print $2 "\t"}' > $ref ; done
