# seqlen.py

## Overview

`seqlen.py` is a Python script designed to generate sequence ID and sequence length pairs from FASTA sequences. It is a direct translation from an `awk` script and is intended to be used in similar contexts where sequence length needs to be calculated quickly and easily.

## Usage

```bash
seqlen.py [fasta_file...]

Arguments

fasta_file: A FASTA file containing one or more sequences. If no file is provided or if the file is specified as '-', the script reads from standard input.

Output

The output is in the format:

<sequence_id> <length>

Example

Given the file Fosmid.fasta containing sequences:

>FFOF1000 3432432 FFOF
ACTG
>FFOF1001 3432433 FFOF
>FFOF1002 3432434 FFOF
ACTG
ACTG

Running the script as follows:

$ python3 seqlen.py Fosmid.fasta

Will produce:

FFOF1000 4
FFOF1001 0
FFOF1002 8

Author

Brandon D. Jordan


### Steps to Implement

1. **Save the Python script** as `seqlen.py`.
2. **Make it executable** (if on a Unix-like system):
   ```bash
   chmod +x seqlen.py

