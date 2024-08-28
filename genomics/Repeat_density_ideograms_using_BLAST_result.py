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

# Define a function to parse arguments
def parse_args():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Generate a chromosome ideogram with BLAST results.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to the input file')
    parser.add_argument('-l', '--chrom_lengths', help="Optional input file with chromosome lengths.")
    parser.add_argument('-w', '--window_size', type=int, help='Window size for density calculation', default=None)
    parser.add_argument('-f', '--filter_file', help='Path to the JSON file containing filter criteria')
    parser.add_argument('-d', '--delimiter', default='\t', help='Delimiter used in the input file (default is tab)')
    parser.add_argument('-c', '--color_break', type=int, default=100000, help='Integer describing how far apart "s. start" values must be to be grouped into a different category. Default is 100 Kb.')

    return parser.parse_args()


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

#    # Check the number of columns
#    num_columns = df.shape[1]
#    print(f"The file has {num_columns} columns after processing the header line.")
#
#    # Check for rows that contain comments
#    rows_with_comments = df[df.apply(lambda row: row.astype(str).str.contains(comment_char).any(), axis=1)]
#
#    if not rows_with_comments.empty:
#        print("The following rows contain comments:")
#        print(rows_with_comments)
#    else:
#        print("No rows with comments found.")

    return df

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
    max_chrom_length = df['s. start'].max()
    return max_chrom_length


def main():
    args = parse_args()

    # Load the filter criteria from the JSON file
    with open(args.filter_file, 'r') as f:  # Replace with your actual file path
        filter_criteria = json.load(f)

#    # Print the loaded filter criteria to verify
#    print("Loaded filter criteria from JSON file:")
#    print(json.dumps(filter_criteria, indent=4))  # Pretty-print the JSON for clarity
#
#    for key, value in filter_criteria.items():
#        print(f"Key: {key}, Value: {value}, Type: {type(value)}")

    df = load_tabular_file(args.input)

    # Load the filter criteria from the JSON file
    with open(args.filter_file, 'r') as f:
        filter_criteria = json.load(f)

    # Subset the DataFrame
    filtered_df = subset_dataframe(df, filter_criteria)

#    print(filtered_df.head())

    # Check the min and max values for the 'alignment length' column
    alignment_length_min = filtered_df['alignment length'].min()
    alignment_length_max = filtered_df['alignment length'].max()

    print(f"Min alignment length: {alignment_length_min}")
    print(f"Max alignment length: {alignment_length_max}")

    # Similarly, you can check other columns as well
    percent_identity_min = filtered_df['% identity'].min()
    percent_identity_max = filtered_df['% identity'].max()

    print(f"Min % identity: {percent_identity_min}")
    print(f"Max % identity: {percent_identity_max}")

def main():
    # Parse the arguments
    args = parse_args()

    if args.filter_file:
        # Load the filter criteria from the JSON file
        with open(args.filter_file, 'r') as f:
            filter_criteria = json.load(f)

    # Load main input file
    try:
        df = load_tabular_file(args.input)
    except KeyError as e:
        print(f"Column not found: {e}")
        raise
    except Exception as e:
        print(f"Error reading the file: {e}")
        raise

    if args.filter_file:
        # Subset the DataFrame
        filtered_df = subset_dataframe(df, filter_criteria)
    else:
        filtered_df = df

    color_break = args.color_break

    # Initialize the 'consolidated_class' column
    filtered_df['consolidated_class'] = ''

    # Initialize variables to track the current group and the last 's. start' value
    current_class = 0
    last_s_start = None
    last_chrom = None

    # First, sort the DataFrame by 'subject acc.ver' and then by 's. start'
    sorted_df = filtered_df.sort_values(by=['subject acc.ver', 's. start'])

    # Iterate over the sorted DataFrame
    for index, row in sorted_df.iterrows():
        current_chrom = row['subject acc.ver']

        # If this is the first row or the chromosome has changed, assign a new class
        if last_chrom is None or current_chrom != last_chrom:
            current_class += 1
            sorted_df.at[index, 'consolidated_class'] = current_class
            last_s_start = row['s. start']
            last_chrom = current_chrom
        else:
            # If within the same chromosome, check the 's. start' difference
            if row['s. start'] - last_s_start > color_break:
                current_class += 1  # Increment the class if the difference exceeds color_break
            sorted_df.at[index, 'consolidated_class'] = current_class
            last_s_start = row['s. start']

#    # QC
#    print(sorted_df[['subject acc.ver','s. start', 'consolidated_class']].head())  # Inspect the first few rows
#    print("Checking for empty 'consolidated_class' values:")
#    print(sorted_df[sorted_df['consolidated_class'] == ''])

    # Convert the consolidated_class column to string to make it easier to work with
    sorted_df['consolidated_class'] = sorted_df['consolidated_class'].astype(str)

    # Extract the unique classes into the selected_classes list
    selected_classes = sorted_df['consolidated_class'].unique()
    print(selected_classes)

    # Initializing chrom_lengths
    chrom_lengths = None

    if args.chrom_lengths:
        chrom_lengths = read_chrom_lengths(args.chrom_lengths)
        max_chrom_length = max(chrom_lengths.values())
        # Unique chromosomes
        chromosomes = chrom_lengths.keys()
    else:
        # If no file is provided, use the original logic for max_chrom_length
        max_chrom_length = original_max_chrom_length_logic(df)
        # Unique chromosomes
        chromosomes = sorted_df['subject acc.ver'].unique()

    # Create a limited palette of 10 colors (or any number you prefer)
    palette = sns.color_palette('rainbow', 10)  # Adjust the number here if you prefer a different size palette

    # Assign colors to each class by cycling through the palette
    class_colors = {cls: palette[i % len(palette)] for i, cls in enumerate(selected_classes)}

    # Create a gradient colormap for each class
    class_colormaps = {}
    for cls, base_color in class_colors.items():
        class_colormaps[cls] = mcolors.LinearSegmentedColormap.from_list(f"{cls}_cmap", ["white", base_color])

    # Create a plot for each chromosome
    fig, axes = plt.subplots(len(chromosomes), 1, figsize=(15, 2.5 * len(chromosomes)), sharex=True)  # Adjusted vertical space

    if len(chromosomes) == 1:
        axes = [axes]

    def calculate_density(chrom_data, chrom_length, window_size):
        density = np.zeros((chrom_length // window_size) + 1)
        for start, end in zip(chrom_data['s. start'], chrom_data['s. end']):
            start_idx = start // window_size
            end_idx = end // window_size
            for i in range(start_idx, end_idx + 1):
                density[i] += 1
        return density

    # Calculate density and max density for each class separately
    density_dict = {}
    max_density_dict = {cls: 0 for cls in selected_classes}

    for chrom in chromosomes:
        chrom_data = sorted_df[sorted_df['subject acc.ver'] == chrom]
        if args.chrom_lengths:
            chrom_length = chrom_lengths.get(chrom)
        else:
            chrom_length = chrom_data['s. start'].max()
        if args.window_size: #
            window_size = args.window_size #
            density_dict[chrom] = {}
            for cls in selected_classes:
                density = calculate_density(chrom_data, chrom_length, window_size)
                density_dict[chrom][cls] = density
                max_density_dict[cls] = max(max_density_dict[cls], max(density))

    for chrom_idx, (ax, chrom) in enumerate(zip(axes, chromosomes)):
        chrom_data = sorted_df[sorted_df['subject acc.ver'] == chrom]
        if args.chrom_lengths:
            chrom_length = chrom_lengths.get(chrom)
        else:
            chrom_length = chrom_data['s. start'].max()
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
                if chrom_idx == 0:
                    ax.text(-0.05 * max_chrom_length, y_position + 0.02, cls, va='center', ha='right', fontsize=10)  # Add class name for the first chromosome
        else:
            # Draw repeat sequences
            for class_idx, cls in enumerate(selected_classes):
                class_data = chrom_data
                y_position = 1.0 - class_idx * 0.12  # Adjusted spacing between tracks
                for _, row in class_data.iterrows():
                    ax.add_patch(mpatches.Rectangle((row['s. start'], y_position), row['s. end'] - row['s. start'], 0.1, color=class_colors[row['consolidated_class']], zorder=1))
                if chrom_idx == 0:
                    ax.text(-0.05 * max_chrom_length, y_position + 0.02, cls, va='center', ha='right', fontsize=10)  # Add class name for the first chromosome
        ax.set_xlim(0, max_chrom_length)
        ax.set_ylim(0.85 - 0.12 * len(selected_classes), 1.2)  # Adjusted ylim to fit tracks properly
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
    output_file = 'Repeat_density_ideogram.svg'
    try:
        fig.subplots_adjust(left=0.1, right=0.85, top=0.95, bottom=0.05)
        plt.savefig(output_file)
        print(f"Plot saved to {output_file}")
    except Exception as e:
        print(f"Error saving the file: {e}")
        raise

    # Show plot (optional, for interactive environments)
    plt.show()

if __name__ == "__main__":
    main()
