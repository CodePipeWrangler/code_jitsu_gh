#!/usr/bin/env python3

import pandas as pd
from collections.abc import MutableMapping
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import argparse
import numpy as np
import matplotlib.colors as mcolors
import json
import itertools

# Define usage function
def print_usage():
    usage_text = """
    Usage: python script_name.py -i [INPUT_FILE] -r [REPEAT_CLASSES] [-o OUTPUT_FILE] [-l CHROM_LENGTHS] [-w WINDOW_SIZE]
                                 [--use_ranges]

    Arguments:
      -i, --input          Path to the input file (required)
      -o, --output         Output filename without extension (optional, default: 'output.svg')
      -r, --repeat_classes Space-separated list of integers or ranges of integers to create classes
                           (e.g., '1-10 20-60') (required)
      -l, --chrom_lengths  Input file with chromosome lengths (optional)
      -w, --window_size    Window size for density calculation. Default is 100kb (optional)
      --use_ranges         Flag to indicate that repeat classes should be treated as ranges (optional)
      -f, --filter_file    Path to the JSON file containing filter criteria (optional)[Temporarily unavailable]
      -d, --delimiter      Delimiter used in the input file (optional, default is tab) [Temporarily unavailable]

    Example:
      python script_name.py -i data.txt -r 1-10 20-30 -o ideogram -l chrom_lengths.txt -w 50000 --use_ranges
    """
    print(usage_text)

# Define a function to parse arguments
def parse_args():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Generate a chromosome ideogram with BLAST results.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to the input file')
    parser.add_argument('-o', '--output', type=str, help='Output filename without extension')
    parser.add_argument(
           '-r', '--repeat_classes', type=str, nargs='+', required=True, help="space-separated list of integers or ranges of integers to create classes (e.g., '1-10 20-60')")
    parser.add_argument('-l', '--chrom_lengths', help="Optional input file with chromosome lengths.")
    parser.add_argument('-w', '--window_size', type=int, help='Window size for density calculation', default=None)
    parser.add_argument('--use_ranges', action='store_true', default=True, help="Set this flag to indicate repeat classes should be treated as ranges. Flag is active by default")
    parser.add_argument('-f', '--filter_file', help='Path to the JSON file containing filter criteria')
    parser.add_argument('-d', '--delimiter', default='\t', help='Delimiter used in the input file (default is tab)') #not currently operational. link with load_tabular_file() in script a few steps downstream
    parser.add_argument('-c', '--color_break', type=int, default=100000, help='Integer describing how far apart "s. start" values must be to be grouped into a       different category. Default is 100 Kb.')

    return parser.parse_args()

# The function below is optimized for BLAST files!!!
def load_tabular_file(file_path, delimiter='\t', comment_char='#', field_identifier='# Fields:'):
    column_names = None
    skip_rows = 0

    with open(file_path, 'r') as file:
        for i, line in enumerate(file):
            if line.startswith(field_identifier):
                # Remove the field identifier and strip whitespace
                header_line = line.replace(field_identifier, '').strip()
                # Split the header by commas to get column names
                column_names = [col.strip() for col in header_line.split(',')]
                skip_rows = i + 1
                break

    if column_names is None:
        raise ValueError("No header line starting with '{}' found in file.".format(field_identifier))

    # Load the data from the file, skipping the lines up to the actual data
    df = pd.read_csv(file_path, delimiter=delimiter, comment=comment_char, names=column_names, skiprows=skip_rows)

    return df

# This should work for ULTRA tab-delimited files. Let's revisit optimizing this later...
def load_ULTRA_file(file_path, num_cols=10, delimiter='\t', comment_char='#'):
    # Read the file with the determined number of columns
    df = pd.read_csv(file_path, delim_whitespace=True, header=None, usecols=range(num_cols), low_memory=False)

    # Assign column names based on the expected structure
    column_names = ["seq_id", "start", "length", "period", "score", "substitutions", "insertions", "deletions", "consensus", "sequence"]
    if num_cols > len(column_names):
        column_names.extend([f'extra_{i}' for i in range(num_cols - len(column_names))])

    df.columns = column_names[:num_cols]

    # Convert relevant columns to numeric types
    df['start'] = pd.to_numeric(df['start'], errors='coerce')
    df['length'] = pd.to_numeric(df['length'], errors='coerce')
    df['period'] = pd.to_numeric(df['period'], errors='coerce')
    df['score'] = pd.to_numeric(df['score'], errors='coerce')
    df['substitutions'] = pd.to_numeric(df['substitutions'], errors='coerce')
    df['insertions'] = pd.to_numeric(df['insertions'], errors='coerce')
    df['deletions'] = pd.to_numeric(df['deletions'], errors='coerce')
    df['seq_id'] = df['seq_id'].astype(str)

    return df

def subset_repeat_classes(df, repeat_classes, use_ranges):
    if use_ranges:
        ranges = []
        for r in repeat_classes:
            if '-' in r:
                start, end = map(int, r.split('-'))
                ranges.append((start, end))
            else:
                # Handle the case where a single integer is provided without a range
                single_value = int(r)
                ranges.append((single_value, single_value))

        df_filtered = pd.concat([df[(df['period'] >= start) & (df['period'] <= end)] for start, end in ranges])
        selected_classes = df_filtered['period'].unique()
    else:
        # Treat repeat_classes as a list of individual integers
        class_list = [int(r) for r in repeat_classes]
        df_filtered = df[df['period'].isin(class_list)]
        selected_classes = df_filtered['period'].unique()

    return df_filtered, selected_classes

def subset_dataframe(df, filter_criteria):
#    # Check the first few rows of the DataFrame
#    print("First few rows of the DataFrame before filtering:")
#    print(df.head())
#
#    # Check the data types of the columns
#    print("\nData types of the columns:")
#    print(df.dtypes)
#
#    # Check the unique values or a summary of 'alignment length'
#    print("\nSummary statistics for 'alignment length' before filtering:")
#    print(df['alignment length'].describe())

    filtered_df = df.copy()

    for column, criteria in filter_criteria.items():
        if column in df.columns:
            if df[column].dtype == 'object':  # Handling string columns
                if isinstance(criteria, str):
                    filtered_df = filtered_df[filtered_df[column].str.contains(criteria)]
            else:  # Handling numeric columns
                if isinstance(criteria, list):
                    if criteria[0] == '>':
                        filtered_df = filtered_df[filtered_df[column] > criteria[1]]
                        #QC
                        print(criteria[0])
                    elif criteria[0] == '<':
                        filtered_df = filtered_df[filtered_df[column] < criteria[1]]
                    elif criteria[0] == 'between':
                        filtered_df = filtered_df[(filtered_df[column] >= criteria[1]) & (filtered_df[column] <= criteria[2])]

    return filtered_df

# Define a function to read chromosome lengths from a file
def read_chrom_lengths(file):
    lengths = {}
    with open(file) as f:
        for line in f:
            chrom, length = line.strip().split()
            lengths[chrom] = int(length)
    return lengths

def original_max_chrom_length_logic(df):
    # Determine the maximum chromosome length for proportional scaling
    max_chrom_length = df['start'].max()
    return max_chrom_length

def main():
    # Parse the arguments
    args = parse_args()

    if not args.input or not args.repeat_classes:
        print_usage()
        exit(1)

    # Load main input file
    try:
        df = load_ULTRA_file(args.input)
    except KeyError as e:
        print(f"Column not found: {e}")
        raise
    except Exception as e:
        print(f"Error reading the file: {e}")
        raise

   # The code immediatley below cannot be run without fitting into the code scheme optimally.
   # if args.filter_file:
   #     # Load the filter criteria from the JSON file
   #     with open(args.filter_file, 'r') as f:
   #         filter_criteria = json.load(f)
   # if args.filter_file:
   #     # Subset the DataFrame
   #     filtered_df = subset_dataframe(df, filter_criteria)
   # else:
   #     filtered_df = df

    # Initializing chrom_lengths
    chrom_lengths = None

    # Map the class/family to the consolidated categories using user-input integers
    df['consolidated_class'] = df['period']

    if args.repeat_classes:
        filtered_df, selected_classes = subset_repeat_classes(df, args.repeat_classes, args.use_ranges)
    #QC
    pass

    # sort the repeat classes so they are printed in numerical order
    selected_classes.sort()

    if args.chrom_lengths:
        chrom_lengths = read_chrom_lengths(args.chrom_lengths)
        max_chrom_length = max(chrom_lengths.values())
        chromosomes = chrom_lengths.keys()
    else:
        max_chrom_length = original_max_chrom_length_logic(filtered_df)
        chromosomes = filtered_df['seq_id'].unique()

    # Assign colors to each class/family using a rainbow palette
    palette = sns.color_palette('rainbow', len(selected_classes))
    class_colors = {cls: palette[i] for i, cls in enumerate(selected_classes)}

    # Create a gradient colormap for each class
    class_colormaps = {}
    for cls, base_color in class_colors.items():
        class_colormaps[cls] = mcolors.LinearSegmentedColormap.from_list(f"{cls}_cmap", ["white", base_color])

    # Number of classes
    num_classes = len(selected_classes)

    # Define a base height per class
    height_per_class = 2.0  # You can adjust this multiplier if needed

    # Create the figure and axes dynamically based on the number of chromosomes and classes
    fig, axes = plt.subplots(len(chromosomes), 1, figsize=(15, 2.5 + num_classes * height_per_class), sharex=True)

    if len(chromosomes) == 1:
        axes = [axes]

    def calculate_density(chrom_data, chrom_length, window_size):
        density = np.zeros((chrom_length // window_size) + 1)
        for start, end in zip(chrom_data['start'], chrom_data['start'] + chrom_data['length']):
            start_idx = start // window_size
            end_idx = end // window_size
            for i in range(start_idx, end_idx + 1):
                density[i] += 1
        return density

    # Calculate density and max density for each class separately
    density_dict = {}
    max_density_dict = {cls: 0 for cls in selected_classes}

    for chrom in chromosomes:
        chrom_data = filtered_df[filtered_df['seq_id'] == chrom]
        if args.chrom_lengths:
            chrom_length = chrom_lengths.get(chrom)
        else:
            chrom_length = chrom_data['start'].max()
        if args.window_size: #
            window_size = args.window_size #
            density_dict[chrom] = {}
            for cls in selected_classes:
                density = calculate_density(chrom_data[chrom_data['consolidated_class'] == cls], chrom_length, window_size)
                density_dict[chrom][cls] = density
                max_density_dict[cls] = max(max_density_dict[cls], max(density))

    for chrom_idx, (ax, chrom) in enumerate(zip(axes, chromosomes)):
        chrom_data = filtered_df[filtered_df['seq_id'] == chrom]
        if args.chrom_lengths:
            chrom_length = chrom_lengths.get(chrom)
        else:
            chrom_length = chrom_data['start'].max()
        # Draw the chromosome ideogram as a horizontal bar
        ax.add_patch(mpatches.Rectangle((0, 1.1), chrom_length, 0.1, color='lightgrey', zorder=0))
        if args.window_size: #
            window_size = args.window_size #
            for class_idx, cls in enumerate(selected_classes):
                density = density_dict[chrom][cls]
                cmap = class_colormaps[cls]
                max_density = max_density_dict[cls]
                norm = mcolors.Normalize(vmin=0, vmax=max_density)
                y_position = 1.0 - class_idx * 0.12  # Adjusted spacing between tracks
                for i, count in enumerate(density):
                    if count > 0:
                        color = cmap(norm(count))
                        ax.add_patch(mpatches.Rectangle((i * window_size, y_position), window_size, 0.1, color=color, zorder=1))
                ax.text(-0.05 * max_chrom_length, y_position + 0.02, cls, va='center', ha='right', fontsize=10)  # Add class name for the first chromosome
        else:
            # Draw repeat sequences
            for class_idx, cls in enumerate(selected_classes):
                class_data = chrom_data[chrom_data['consolidated_class'] == cls]
                y_position = 1.0 - class_idx * 0.12  # Adjusted spacing between tracks
                for _, row in class_data.iterrows():
                    ax.add_patch(mpatches.Rectangle((row['start'], y_position), row['start'] + row['length'] - row['start'], 0.1, color=class_colors[row['consolidated_class']], zorder=1))
                ax.text(-0.05 * max_chrom_length, y_position + 0.02, cls, va='center', ha='right', fontsize=10)  # Add class name for the first chromosome
        ax.set_xlim(0, max_chrom_length)
        ax.set_ylim(0.85 - num_classes * 0.12, 1.2)  # Adjusted ylim to fit tracks properly
        ax.set_yticks([])
        ax.set_title(f'Chromosome: {chrom} (Size: {chrom_length} bp)', fontsize=10)

    # Create a legend for repeat classes and position it to avoid overlap with the density legend
    handles = [mpatches.Patch(color=class_colors[cls], label=cls) for cls in selected_classes]
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.05, 1))

    # Add colorbars for density of each repeat class
    if args.window_size: #
        plt.subplots_adjust(right=0.85)  # Adjust right to make space for colorbars
        for idx, cls in enumerate(selected_classes):
            max_density = max_density_dict[cls]
            norm = mcolors.Normalize(vmin=0, vmax=max_density)
            cax = fig.add_axes([0.91, 0.1 + idx * 0.1, 0.02, 0.08])
            sm = plt.cm.ScalarMappable(cmap=class_colormaps[cls], norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, cax=cax)
            cbar.set_label(f'Density of {cls}')

    # Save the plot to a SVG file
    if args.output:
        output_file = f"{args.output}.svg"
    else:
        output_file = 'Repeat_density_ideogram.svg'

    try:
        fig.subplots_adjust(left=0.1, right=0.85, top=0.95, bottom=0.05)
        plt.savefig(output_file)
        print(f"Plot saved to {output_file}")
    except Exception as e:
        print(f"Error saving the file: {e}")
        raise

    # Show plot (optional, for interactive environments)
#    plt.show()

if __name__ == "__main__":
    main()


