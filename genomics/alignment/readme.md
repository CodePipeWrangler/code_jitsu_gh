# Sequence Alignment Script

This script performs multiple sequence alignment using Clustal Omega on sequences provided in a FASTA file. The aligned sequences are saved in a new FASTA file, and the alignment can be further processed or converted to other formats.

## Requirements

- Python 3.x
- [Biopython](https://biopython.org/) library
- Clustal Omega executable

### Installation

To install the required Python library, use:

```bash
pip install biopython

Ensure that Clustal Omega is installed and accessible in your system's PATH. You can download it from the Clustal Omega website.

Usage
python align_fasta.py <input_fasta> <output_fasta> [clustalo_path]

Arguments
<input_fasta>: Path to the input FASTA file containing sequences to be aligned.
<output_fasta>: Path to the output FASTA file where the aligned sequences will be saved.
[clustalo_path]: Optional path to the Clustal Omega executable. If omitted, the script assumes clustalo is available in your system's PATH.
Example
Align sequences from input_sequences.fasta and save the result to aligned_sequences.fasta:

python align_fasta.py input_sequences.fasta aligned_sequences.fasta

If Clustal Omega is installed in a non-standard location, specify the path:

python align_fasta.py input_sequences.fasta aligned_sequences.fasta /custom/path/to/clustalo

Output

The script will output the aligned sequences in the specified output FASTA file. It will also print the alignment to the console.

Optionally, the alignment can be saved in other formats such as PHYLIP:

AlignIO.write(alignment, "aligned_sequences.phy", "phylip")

Troubleshooting

If you encounter a FileNotFoundError for Clustal Omega, ensure that the executable is correctly installed and accessible via the provided path or the system's PATH.

License

This project is licensed under the MIT License. See the LICENSE file for details.

Contributing

Contributions are welcome! Please fork this repository and submit a pull request.

Contact

