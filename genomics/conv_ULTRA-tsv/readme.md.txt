# ultra2tsv.sh

This script is designed to convert the JSON output of the ULTRA tool (https://github.com/TravisWheelerLab/ULTRA) into TSV (Tab-Separated Values) format. ULTRA is a software tool for identifying and characterizing repetitive elements in DNA sequences, particularly useful for analyzing centromeric repeats and transposable elements.

## Features

- **Converts JSON to TSV:** The script processes ULTRA's JSON output files and converts them to a more accessible TSV format.
- **Targeted Data Extraction:** It captures key fields from the JSON, including `SequenceName`, `Start`, `Length`, `Period`, `Score`, `Substitutions`, `Insertions`, `Deletions`, `Consensus`, and `Sequence`.
- **Usage Instructions:** The script includes a usage function that provides instructions when run without parameters or with the `--help` option.

## Under Construction

Currently, the script does not capture the Sub-repeat data. This functionality is under development. However, most of the useful data for identifying centromeric repeats or transposable elements is already embedded in the Repeat data alone, making the script highly useful even in its current form.

## Usage

To run the script, use the following command:

```bash
./ultra2tsv.sh <json_file>

Replace `<json_file>` with the path to the JSON file you wish to convert. The script will generate a TSV file with the same name as the JSON file but with a `.tsv` extension.

If you need help or want to see usage instructions, you can run:

```bash
./ultra2tsv.sh --help
```

## Requirements

- **Bash Shell:** This script is intended to be run in a Unix-like environment with a Bash shell.
- **Perl:** The script uses Perl for processing JSON data.
- **awk:** Utilized for extracting specific fields from the JSON data.

## Future Work

- **Sub-repeat Data Capture:** A future version of this script will include functionality to capture and convert Sub-repeat data from the JSON files.

## Contact

For any questions or suggestions, please contact Brandon Jordan at brandon.jordan@usda.gov.
