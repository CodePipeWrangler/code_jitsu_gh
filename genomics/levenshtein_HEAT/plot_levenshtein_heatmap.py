#!/usr/bin/env python3

import argparse
import re
import Levenshtein as pylev
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO

def remove_pattern(seq_id, pattern):
    """Remove a specified pattern from the sequence ID."""
    return re.sub(pattern, '', seq_id)

def truncate_id(seq_id, length=15):
    """Truncate the sequence ID to a specified length."""
    return seq_id[:length]

def process_ids(seqs, pattern=None):
    """Process sequence IDs by removing a pattern and then truncating."""
    if pattern:
        seqs = {remove_pattern(seq_id, pattern): seq for seq_id, seq in seqs.items()}
    seqs = {truncate_id(seq_id): seq for seq_id, seq in seqs.items()}
    return seqs

def calculate_distances(seqs, label_type):
    """Calculate Levenshtein distance between sequences."""
    seq_ids = list(seqs.keys())
    chain_r = {}

    for i in range(len(seq_ids)):
        for j in range(i, len(seq_ids)):
            if i == j:
                continue
            id1, id2 = seq_ids[i], seq_ids[j]
            seq1, seq2 = seqs[id1], seqs[id2]
            distance = pylev.distance(seq1, seq2)
            ratio = pylev.ratio(seq1, seq2)

            if id1 not in chain_r:
                chain_r[id1] = {}
            chain_r[id1][id2] = ratio

            if id2 not in chain_r:
                chain_r[id2] = {}
            chain_r[id2][id1] = ratio

    return pd.DataFrame(chain_r).sort_index().sort_index(axis=1)

def usage():
    """Prints usage instructions."""
    usage_message = """
Usage:
    plot_levenshtein_heatmap.py [-h] fasta_file [-l {id,sequence}] [-p pattern]

Arguments:
    fasta_file          Input FASTA file containing sequences.

Optional arguments:
    -l, --labels        Choose whether to display sequence ID or sequence on the vertical axis of the plot
                        (default: id).
    -p, --pattern       Pattern to remove from sequence IDs before processing. Enclose in quotes.

Description:
    This script calculates Levenshtein distances between sequences in a FASTA file and generates a heatmap.
    It can be used in genetics and other applications where string distance measurements are needed.
    """
    print(usage_message)

def main():
    parser = argparse.ArgumentParser(description="Calculate Levenshtein distances from a FASTA file and plot a heatmap.", add_help=False)
    parser.add_argument("fasta_file", nargs="?", help="Input FASTA file containing sequences")
    parser.add_argument("-l", "--labels", choices=["id", "sequence"], default="id",
                        help="Choose whether to display sequence ID or sequence on the vertical axis of the plot (default: id)")
    parser.add_argument("-p", "--pattern", help="Pattern to remove from sequence IDs before processing (enclose in quotes)")

    args = parser.parse_args()

    if not args.fasta_file:
        usage()
        return

    seqs = {record.id: str(record.seq) for record in SeqIO.parse(args.fasta_file, "fasta")}

    if args.pattern:
        seqs = process_ids(seqs, pattern=args.pattern)
    else:
        seqs = process_ids(seqs)

    df_ratios = calculate_distances(seqs, args.labels)

    sns.heatmap(df_ratios)
    plt.show()

if __name__ == "__main__":
    main()

