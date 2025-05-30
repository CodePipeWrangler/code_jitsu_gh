#!/usr/bin/env python3

# Tandem Repeat Data Statistics
# Calculate and visualize GC content

# Import libraries
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
from sklearn.preprocessing import StandardScaler

# Sequence-based features: GC content
def gc_content(seq_column, seq_column_name):
    """
    Calculate GC content from a specified column in a dataframe.
    Arguments:
    - seq_column: The column of sequences in the dataframe (string values).
    - seq_column_name: The name of the sequence column.
    Returns:
    - A pandas series with the GC content for each sequence.
    """
    return seq_column.apply(lambda seq: (seq.upper().count("G") 
                                         + seq.upper().count("C")) / len(seq))


# Read the file and return the dataframe
def read_file(filename):
    """
    Reads the ULTRA file and returns a dataframe.
    """
    return pd.read_csv(filename, sep="\t", names=['chr', 'pos', 'len', 
                        'period', 'score', 'sub', 'ins', 'del', 'cons', 
                        'seq'], usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

# Normalize or scale the feature matrix
def normalize_features(df, feature_columns):
    """
    Normalizes the feature columns in the dataframe.
    Arguments:
    - df: DataFrame containing the data.
    - feature_columns: List of columns to normalize.
    Returns:
    - Scaled numpy array.
    """
    X_numeric = df[feature_columns].values
    scaler = StandardScaler()
    return scaler.fit_transform(X_numeric)

# Create bins from positions
def create_bins(df, bin_size=1_000_000, position_column="pos"):
    """
    Creates bins in the dataframe from the position column.
    Arguments:
    - df: DataFrame containing the data.
    - bin_size: Size of each bin (default is 1,000,000).
    - position_column: The column to be used for binning.
    """
    df["bin"] = (df[position_column] // bin_size) * bin_size  # start of bin

# Group by bin and compute mean GC content
def group_by_bin(df, bin_column="bin", gc_column="gc_content"):
    """
    Groups the dataframe by the bin and calculates mean GC content.
    Arguments:
    - df: DataFrame containing the data.
    - bin_column: The column used for binning.
    - gc_column: The column containing GC content.
    Returns:
    - A dataframe with mean GC content per bin.
    """
    return df.groupby(bin_column)[gc_column].mean().reset_index()

# Optional: Count number of repeats per bin
def count_repeats_per_bin(df, gc_by_bin, bin_column="bin", 
                          gc_column="gc_content"):
    """
    Counts the number of repeats per bin.
    Arguments:
    - df: DataFrame containing the data.
    - gc_by_bin: DataFrame with GC content per bin.
    - bin_column: The column used for binning.
    - gc_column: The column containing GC content.
    Returns:
    - A dataframe with the repeat count for each bin.
    """
    gc_by_bin["repeat_count"] = df.groupby(bin_column)[gc_column].count().values
    return gc_by_bin

# Plot GC content by bin and save the figure
def plot_gc_content_by_bin(gc_by_bin, bin_size, i, output_filename="output_plot.png"):
    """
    Plots the GC content by bin and saves the figure.
    Arguments:
    - gc_by_bin: DataFrame with GC content per bin.
    - bin_size: The size of each bin.
    - i: The index of the plot (to make unique file names).
    - output_filename: The name of the output plot file 
    (default is "output_plot.png").
    """
    plt.figure(figsize=(12, 5))
    plt.plot(gc_by_bin["bin"], gc_by_bin["gc_content"], marker="o")
    plt.title(f"Mean GC Content per {bin_size:,} bp Bin")
    # scale X-axis manually
    plt.xticks(ticks=gc_by_bin["bin"][::10], 
               labels=(gc_by_bin["bin"][::10] // 1_000_000))
    plt.xlabel("Genomic Position (Mb)")
    plt.ylabel("Mean GC Content")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    plt.cla()

# Main function to process files
def process_files_plot_GC_by_bin(num_files, filename_pattern, seq_column_name, feature_columns, bin_size=1_000_000):
    """
    Main function to process multiple files with the given parameters.
    Arguments:
    - num_files: Number of files to process.
    - filename_pattern: Pattern of filenames (with '#' for file number).
    - seq_column_name: The name of the column containing sequences.
    - feature_columns: List of feature columns for normalization.
    - bin_size: The size of bins for GC content calculation 
    (default is 1,000,000).
    """
    for i in range(1, num_files + 1):
        filename = filename_pattern.replace("#", str(i))
        print(f"Processing file: {filename}")

        # Step 1: Read the file
        df = read_file(filename)

        # Step 2: Normalize features
        X_scaled = normalize_features(df, feature_columns)

        # Step 3: Calculate GC content
        df["gc_content"] = gc_content(df[seq_column_name], seq_column_name)

        # Step 4: Add GC content to the feature matrix
        X_full = np.hstack([X_scaled, df[["gc_content"]].values])

        # Step 5: Create bins
        create_bins(df, bin_size)

        # Step 6: Group by bin and calculate mean GC content
        gc_by_bin = group_by_bin(df)

        # Step 7: Optional: Count repeats per bin
        gc_by_bin = count_repeats_per_bin(df, gc_by_bin)

        # Step 8: Plot and save
        plot_gc_content_by_bin(gc_by_bin, bin_size, i, 
                               f"output_{i}_gc_content_by_bin.png")