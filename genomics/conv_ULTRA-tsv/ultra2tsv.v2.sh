#!/bin/bash

# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=05:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=2   # 1 processor core(s) per node X 2 threads per core
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="ultra2tsv"
#SBATCH --mail-user=brandon.jordan@usda.gov   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

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
