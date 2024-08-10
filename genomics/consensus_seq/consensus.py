#!/usr/bin/env python

import sys
from Bio import AlignIO
from Bio.Align import AlignInfo

def get_consensus_sequence(alignment_file, consensus_threshold):
    alignment = AlignIO.read(alignment_file, 'fasta')
    summary_align = AlignInfo.SummaryInfo(alignment)
    return summary_align.dumb_consensus(float(consensus_threshold))

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <alignment_file> <consensus_threshold>")
        sys.exit(1)

    alignment_file = sys.argv[1]
    consensus_threshold = sys.argv[2]

    consensus_sequence = get_consensus_sequence(alignment_file, consensus_threshold)
    print(consensus_sequence)

if __name__ == "__main__":
    main()

