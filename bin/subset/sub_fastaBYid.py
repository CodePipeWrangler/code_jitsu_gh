#!/usr/bin/env python3

import sys
from Bio import SeqIO

def usage():
    print("""
Usage: python subset_fasta.py <input_fasta> <id_list_file> <output_fasta>

Arguments:
    <input_fasta>   : Path to the input FASTA file.
    <id_list_file>  : Path to the file containing sequence IDs to keep, one per line.
    <output_fasta>  : Path to the output FASTA file that will contain the subset sequences.

Example:
    python subset_fasta.py input.fasta id_list.txt subset.fasta
""")
    sys.exit(1)

def subset_fasta(input_fasta, id_list_file, output_fasta):
    # Load the list of IDs to keep
    with open(id_list_file) as f:
        ids_to_keep = set(line.strip() for line in f)

    # Open the output file for writing
    with open(output_fasta, "w") as output_handle:
        # Iterate through the FASTA file and write matching sequences to the output file
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id in ids_to_keep:
                SeqIO.write(record, output_handle, "fasta")

    print(f"Subset FASTA file created: {output_fasta}")

if __name__ == "__main__":
    # Check if the correct number of arguments are provided
    if len(sys.argv) != 4:
        usage()

    # Get arguments from the command line
    input_fasta = sys.argv[1]
    id_list_file = sys.argv[2]
    output_fasta = sys.argv[3]

    # Run the subset function
    subset_fasta(input_fasta, id_list_file, output_fasta)

