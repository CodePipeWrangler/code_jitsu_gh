#!/usr/bin/python3

import sys
from Bio import SeqIO

def usage():
    print("Usage: python analyze_consensus.py <sequence> or <fasta_file>")
    print("Examples:")
    print("  python analyze_consensus.py AAAACGAAGCAACAnGCATCTTCCCCTCAACTCTAACCTAAGATACCATTTAATTACTTG...")
    print("  python analyze_consensus.py consensus_sequences.fasta")
    sys.exit(1)

def analyze_sequence(sequence, cluster_threshold=5):
    # Calculate the total length of the sequence
    seq_length = len(sequence)

    # Count the number of 'N' characters
    n_count = sequence.upper().count('N')

    # Calculate the proportion of 'N' characters
    n_proportion = n_count / seq_length if seq_length > 0 else 0

    # Find positions of 'N' characters (1-based index)
    n_positions = [i + 1 for i, base in enumerate(sequence.upper()) if base == 'N']

    # Identify clustered 'N' characters
    clusters = []
    current_cluster = [n_positions[0]] if n_positions else []

    for i in range(1, len(n_positions)):
        if n_positions[i] - n_positions[i-1] <= cluster_threshold:
            current_cluster.append(n_positions[i])
        else:
            if len(current_cluster) > 1:
                clusters.append(current_cluster)
            current_cluster = [n_positions[i]]

    if len(current_cluster) > 1:  # Append the last cluster if it has more than one element
        clusters.append(current_cluster)

    # Print the results
    print(f"Sequence Length: {seq_length} bp")
    print(f"Number of ambiguous calls ('N'): {n_count}")
    print(f"Proportion of ambiguous calls: {n_proportion:.4f} ({n_proportion * 100:.2f}%)")
    print(f"Positions of ambiguous calls: {n_positions}")

    if clusters:
        print(f"Ambiguous calls are clustered (within {cluster_threshold} bases).")
        for cluster in clusters:
            print(f"Clustered ambiguous calls at positions: {cluster}")
    else:
        print("Ambigous calls are not clustered.")

def analyze_fasta(fasta_file, cluster_threshold=5):
    for record in SeqIO.parse(fasta_file, "fasta"):
        print(f"\nAnalyzing sequence {record.id}:")
        analyze_sequence(str(record.seq), cluster_threshold)

def main():
    if len(sys.argv) < 2:
        usage()

    input_data = sys.argv[1]
    cluster_threshold = 5  # Default cluster threshold

    if input_data.endswith(".fasta") or input_data.endswith(".fa"):
        # Treat input as a FASTA file
        analyze_fasta(input_data, cluster_threshold)
    else:
        # Treat input as a direct sequence
        print("Analyzing provided sequence:")
        analyze_sequence(input_data, cluster_threshold)

if __name__ == "__main__":
    main()

