#!/usr/bin/env python3

import sys
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
import os

def usage():
    print("Usage: python align_fasta.py <input_fasta> <output_fasta> [clustalo_path]")
    print("Arguments:")
    print("  <input_fasta>   : Path to the input FASTA file containing sequences to be aligned.")
    print("  <output_fasta>  : Path to the output FASTA file where the aligned sequences will be saved.")
    print("  [clustalo_path] : Optional path to the Clustal Omega executable. Defaults to 'clustalo'.")
    sys.exit(1)

def align_sequences(input_fasta, output_fasta, clustalo_path="clustalo"):
    # Define the command line for Clustal Omega
    clustalomega_cline = ClustalOmegaCommandline(
        cmd=clustalo_path,
        infile=input_fasta,
        outfile=output_fasta,
        verbose=True,
        auto=True
    )

    # Run the Clustal Omega alignment
    stdout, stderr = clustalomega_cline()

    # Parse the output aligned sequences
    alignment = AlignIO.read(output_fasta, "fasta")

    return alignment

def main():
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        usage()

    # Get arguments
    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    clustalo_path = sys.argv[3] if len(sys.argv) == 4 else "clustalo"

    # Perform the alignment
    alignment = align_sequences(input_fasta, output_fasta, clustalo_path)

    # Display the alignment
    print(alignment)

    # Optionally, you can save the alignment to another format or process it further
    AlignIO.write(alignment, "aligned_sequences.phy", "clustal")

if __name__ == "__main__":
    main()

