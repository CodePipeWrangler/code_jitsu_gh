#!/usr/bin/env python3

import os
import sys
from Bio import AlignIO
from Bio.SeqIO import write
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.Align.Applications import ClustalOmegaCommandline

def usage():
    print("Usage: python optimize_workflow.py <input_fasta> <output_fasta> [clustalo_path] [max_x] [trim_threshold]")
    sys.exit(1)

def align_sequences(input_fasta, output_fasta, clustalo_path="clustalo"):
    clustalomega_cline = ClustalOmegaCommandline(
        cmd=clustalo_path,
        infile=input_fasta,
        outfile=output_fasta,
        verbose=True,
        auto=True,
        force=True # Force overwrite
    )
    stdout, stderr = clustalomega_cline()
    alignment = AlignIO.read(output_fasta, "fasta")
    return alignment

def trim_alignment(alignment, threshold=0.5):
    trimmed = []
    for record in alignment:
        gaps = record.seq.count('-')
        if gaps / len(record.seq) < threshold:
            trimmed.append(record)
    return MultipleSeqAlignment(trimmed)

def generate_consensus(alignment):
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus(threshold=0.7, ambiguous='X')
    return str(consensus)

def main(input_fasta, output_fasta, max_x=5, clustalo_path="clustalo", trim_threshold=0.5):
    iteration = 0
    while True:
        iteration += 1
        print(f"Iteration {iteration}: Aligning sequences...")
        alignment = align_sequences(input_fasta, output_fasta, clustalo_path)

        print("Trimming sequences...")
        trimmed_alignment = trim_alignment(alignment, threshold=trim_threshold)

        print("Generating consensus...")
        consensus_sequence = generate_consensus(trimmed_alignment)
        print(f"Consensus: {consensus_sequence}")

        if consensus_sequence.count('X') <= max_x:
            print("Consensus sequence is satisfactory.")
            break

        print("Consensus contains too many 'X' characters, re-aligning with trimmed sequences.")

        trimmed_fasta = "trimmed_sequences.fasta"
        write(trimmed_alignment, trimmed_fasta, "fasta")
        input_fasta = trimmed_fasta

if __name__ == "__main__":
    if len(sys.argv) < 3:
        usage()

    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    clustalo_path = sys.argv[3] if len(sys.argv) > 3 else "clustalo"
    max_x = int(sys.argv[4]) if len(sys.argv) > 4 else 5
    trim_threshold = float(sys.argv[5]) if len(sys.argv) > 5 else 0.5

    main(input_fasta, output_fasta, max_x=max_x, clustalo_path=clustalo_path, trim_threshold=trim_threshold)

