import sys
import pandas as pd
import matplotlib.pyplot as plt

def load_data(file_path):
    """
    Load the DNA repeat data from a file into a pandas DataFrame.

    Parameters:
    file_path (str): The path to the data file.

    Returns:
    pd.DataFrame: The loaded data as a DataFrame.
    """
    df = pd.read_csv(file_path, sep='\t', header=None)
    return df

def plot_histogram(df, column_index, highlight_value=158, bins=50, title='Histogram of Repeat Sizes', xlabel='Repeat Size', ylabel='Frequency', repeat_range=None, frequency_range=None):
    """
    Plot a dotted histogram of the repeat sizes with color coding for a specific value.

    Parameters:
    df (pd.DataFrame): The DataFrame containing the DNA repeat data.
    column_index (int): The index of the column containing the repeat sizes.
    highlight_value (int): The value to be highlighted in a different color. Default is 158.
    bins (int): The number of bins for the histogram. Default is 50.
    title (str): The title of the plot. Default is 'Histogram of Repeat Sizes'.
    xlabel (str): The label for the x-axis. Default is 'Repeat Size'.
    ylabel (str): The label for the y-axis. Default is 'Frequency'.
    repeat_range (tuple): A tuple specifying the (min, max) range of repeat sizes to display. Default is None (no limit).
    frequency_range (tuple): A tuple specifying the (min, max) range of frequencies to display on the y-axis. Default is None (no limit).
    """
    # Separate the data into highlighted and non-highlighted parts
    data = df.iloc[:, column_index]

    # Apply the range filter if specified
    if repeat_range is not None:
        data = data[(data >= repeat_range[0]) & (data <= repeat_range[1])]

    highlight_data = data[data == highlight_value]
    other_data = data[data != highlight_value]

    # Plot the histogram for non-highlighted data as dotted lines
    plt.figure(figsize=(10, 6))
    plt.hist(other_data, bins=bins, histtype='step', linestyle=':', edgecolor='black', linewidth=1.5, alpha=0.7, label='Other Sizes', range=repeat_range)

    # Overlay the histogram for highlighted data as dotted lines with different color
    if not highlight_data.empty:
        plt.hist(highlight_data, bins=bins, histtype='step', linestyle=':', edgecolor='red', linewidth=1.5, alpha=0.7, label=f'Size {highlight_value}', range=repeat_range)

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # Apply the frequency range if specified
    if frequency_range is not None:
        plt.ylim(frequency_range)

    plt.legend()
    plt.grid(True)
    plt.show()

def usage():
    print("Usage: python histogram_plot.py <file_path> [options]")
    print("Options:")
    print("  -h, --help       Show this help message and exit")
    print("  --bins=BINS      Number of bins for the histogram (default=50)")
    print("  --highlight=VAL  Highlight this specific repeat size (default=158)")
    print("  --range=MIN,MAX  Range of repeat sizes to include (e.g., 100,500)")
    print("  --freq_range=MIN,MAX  Range of frequencies to display on the y-axis (e.g., 0,50)")
    print("\nExample:")
    print("  python histogram_plot.py data.txt --bins=30 --highlight=150 --range=100,200 --freq_range=0,10")

def main(args):
    if len(args) < 2 or '-h' in args or '--help' in args:
        usage()
        sys.exit()

    file_path = args[1]
    # Assume rest of the code processes arguments and plots histogram

if __name__ == "__main__":
    main(sys.argv)

