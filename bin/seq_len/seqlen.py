#!/usr/bin/env python3

import sys

def calculate_sequence_lengths(file_handle):
    sequence_id = None
    sequence_length = 0

    for line in file_handle:
        line = line.strip()
        if line.startswith('>'):
            if sequence_id is not None:
                print(f"{sequence_id}\t{sequence_length}")
            sequence_id = line[1:].split()[0]
            sequence_length = 0
        else:
            sequence_length += len(line)

    # Print results for the last sequence
    if sequence_id is not None:
        print(f"{sequence_id}\t{sequence_length}")

def usage():
    print(f"Usage: {sys.argv[0]} [fasta_file...]")
    print("Generate sequence ID & sequence length from FASTA sequence(s).")
    print("\nIf no input file is specified or if the file is '-', reads from standard input.")
    sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        usage()

    for fasta_file in sys.argv[1:]:
        if fasta_file == '-':
            calculate_sequence_lengths(sys.stdin)
        else:
            with open(fasta_file, 'r') as file_handle:
                calculate_sequence_lengths(file_handle)

