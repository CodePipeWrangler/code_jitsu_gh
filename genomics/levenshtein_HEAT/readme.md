Levenshtein Distance Heatmap Generator

This script calculates Levenshtein distances between sequences in a FASTA file and generates a heatmap. The script is designed to handle both sequence IDs and the actual sequences for labeling the heatmap. It also includes an option to remove a specified pattern from sequence IDs before processing, making it flexible for different data formats and naming conventions.

Features

Levenshtein Distance Calculation: Computes both the distance and ratio between pairs of sequences.
Heatmap Visualization: Generates a heatmap to visually represent the similarity or dissimilarity between sequences.
Flexible Labeling: Allows users to choose whether to label the heatmap by sequence IDs or the sequences themselves.
Pattern Removal: Users can remove specific patterns from sequence IDs before processing to ensure distinct and meaningful labels.

Usage

python plot_levenshtein_heatmap.py [-h] fasta_file [-l {id,sequence}] [-p pattern]

Arguments

fasta_file: Input FASTA file containing sequences.

Optional Arguments

-l, --labels: Choose whether to display sequence ID or sequence on the vertical axis of the plot (default: id).
-p, --pattern: Pattern to remove from sequence IDs before processing. Enclose in quotes (e.g., 'GENUS.SPECIES.genome.').

Example

To generate a heatmap using sequence IDs as labels and remove a specific pattern from the IDs:

python plot_levenshtein_heatmap.py sequences.fasta --labels id --pattern 'bauva.BV-YZ2020.gnm2.'

This command will remove the pattern 'bauva.BV-YZ2020.gnm2.' from each sequence ID, truncate the remaining part of the ID to 15 characters, and then generate the heatmap.

Applications

While this script was developed for genetic analysis, particularly in comparative genomics, it is versatile enough for any field that requires measuring string distances. Examples include:

Bioinformatics: Analyzing sequence similarity in genetic data.
Text Processing: Comparing document similarity based on string distances.
Data Deduplication: Identifying similar entries in datasets.

Installation

Ensure you have the required dependencies installed:

pip install biopython python-Levenshtein seaborn

License

This project is licensed under the MIT License - see the LICENSE file for details.

Contributing

If you have suggestions for improvements or find bugs, feel free to open an issue or submit a pull request.

Author

Dr. Brandon D. Jordan
