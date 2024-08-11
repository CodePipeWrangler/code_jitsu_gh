# DNA Repeat Data Histogram Plotter

## Overview
This Python script generates histogram plots for DNA repeat data, ideal for publication-quality visuals. It features options to highlight specific repeat sizes, adjust the number of bins, and limit the display range for both repeat sizes and frequencies.

## Requirements
- Python 3
- pandas
- matplotlib

## Usage
Run the script from the command line with the path to your data file and optional parameters to customize the histogram:

### Options
- `-h, --help`: Show help message and exit.
- `--bins=BINS`: Set the number of bins for the histogram (default is 50).
- `--highlight=VAL`: Highlight entries with this repeat size (default is 158).
- `--range=MIN,MAX`: Specify the range of repeat sizes to include (e.g., 100,500).
- `--freq_range=MIN,MAX`: Set the y-axis range for frequencies (e.g., 0,50).

### Example


## Contributing
Feel free to fork this repository and submit pull requests or open an issue for any bugs or enhancements.

## License
Distributed under the MIT License. See `LICENSE` for more information.


Author

Dr. Brandon D. Jordan
