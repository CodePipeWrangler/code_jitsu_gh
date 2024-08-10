import pandas as pd
from collections.abc import MutableMapping
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import argparse
import numpy as np
import matplotlib.colors as mcolors

# Define a function to map repeat classes to consolidated categories
def map_repeat_class(class_family):
    if class_family.startswith('DNA/'):
        return 'DNA transposon'
    elif class_family.startswith('LTR/'):
        return 'LTR elements'
    elif class_family == 'Simple_repeat':
        return 'Simple_repeat'
    elif class_family == 'Low_complexity':
        return 'Low_complexity'
    elif class_family in ['snRNA', 'tRNA', 'rRNA', 'scRNA']:
        return 'Small RNA'
    elif class_family.startswith('RC/'):
        return 'Rolling-circles'
    elif class_family.startswith('LINE/'):
        return 'LINE'
    elif class_family.startswith('SINE/'):
        return 'SINE'
    elif class_family == 'Unknown':
        return 'Unknown'
    else:
        return 'Other'  # If there's any other category not mentioned

# Set up argument parser
parser = argparse.ArgumentParser(description='Generate a chromosome ideogram with repeat sequences.')
parser.add_argument('-i', '--input', type=str, required=True, help='Path to the input file')
parser.add_argument(
    '-r', '--repeat_classes', type=str, default='1,2,3,4,5,6,7,8,9',
    help=(
        'Comma-separated list of repeat classes to plot (1-9):\n'
        '1: DNA transposon\n'
        '2: LTR elements\n'
        '3: Simple_repeat\n'
        '4: Low_complexity\n'
        '5: Small RNA\n'
        '6: Rolling-circles\n'
        '7: LINE\n'
        '8: SINE\n'
        '9: Unknown'
    )
)
parser.add_argument('-w', '--window_size', type=int, help='Window size for density calculation', default=None)
args = parser.parse_args()

# Read the input file
file_path = args.input
try:
    # Dynamically determine the number of columns
    with open(file_path, 'r') as f:
        first_line = f.readline()
        num_cols = len(first_line.split())

    # Read the file with the determined number of columns
    df = pd.read_csv(file_path, delim_whitespace=True, header=None, usecols=range(num_cols), low_memory=False)

    # Assign column names based on the expected structure, add dummy names for extra columns
    column_names = ["score", "perc_div", "perc_del", "perc_ins", "sequence", "query_begin", "query_end",
                    "query_left", "strand", "repeat", "class_family", "repeat_begin", "repeat_end"]
    if num_cols > len(column_names):
        column_names.extend([f'extra_{i}' for i in range(num_cols - len(column_names))])
    df.columns = column_names[:num_cols]

    # Remove the first two header rows
    df = df.drop([0, 1]).reset_index(drop=True)

    # Print the column names to verify
    print("Column names after reading the file:", df.columns.tolist())

    # Convert relevant columns to numeric types
    df['query_begin'] = pd.to_numeric(df['query_begin'], errors='coerce')
    df['query_end'] = pd.to_numeric(df['query_end'], errors='coerce')
    df['query_left'] = pd.to_numeric(df['query_left'].str.replace(r'\(', '', regex=True).str.replace(r'\)', '', regex=True), errors='coerce')
    df['repeat_begin'] = pd.to_numeric(df['repeat_begin'], errors='coerce')
    df['repeat_end'] = pd.to_numeric(df['repeat_end'], errors='coerce')

    # Check if 'repeat_left' exists and process it if it does
    if 'repeat_left' in df.columns:
        df['repeat_left'] = pd.to_numeric(df['repeat_left'].str.replace(r'\(', '', regex=True).str.replace(r'\)', '', regex=True), errors='coerce')
    else:
        print("'repeat_left' column not found, skipping this column.")

except KeyError as e:
    print(f"Column not found: {e}")
    raise
except Exception as e:
    print(f"Error reading the file: {e}")
    raise

# Print the DataFrame to verify the data
print(df.head())

# Map the class/family to the consolidated categories
df['consolidated_class'] = df['class_family'].apply(map_repeat_class)

# Define the mapping from numeric options to consolidated categories
class_mapping = {
    '1': 'DNA transposon',
    '2': 'LTR elements',
    '3': 'Simple_repeat',
    '4': 'Low_complexity',
    '5': 'Small RNA',
    '6': 'Rolling-circles',
    '7': 'LINE',
    '8': 'SINE',
    '9': 'Unknown'
}

# Parse the repeat classes to plot
selected_classes = [class_mapping[c] for c in args.repeat_classes.split(',')]

# Ensure the classes are in the specified order
class_order = ['DNA transposon', 'LTR elements', 'Simple_repeat', 'Low_complexity', 'Small RNA', 'Rolling-circles', 'LINE', 'SINE', 'Unknown']
selected_classes = [cls for cls in class_order if cls in selected_classes]

# Filter the DataFrame to include only the selected repeat classes
df = df[df['consolidated_class'].isin(selected_classes)]

# Unique chromosomes
chromosomes = df['sequence'].unique()

# Assign colors to each class/family using a rainbow palette
palette = sns.color_palette('rainbow', len(selected_classes))
class_colors = {cls: palette[i] for i, cls in enumerate(selected_classes)}

# Create a gradient colormap for each class
class_colormaps = {}
for cls, base_color in class_colors.items():
    class_colormaps[cls] = mcolors.LinearSegmentedColormap.from_list(f"{cls}_cmap", ["white", base_color])

# Determine the maximum chromosome length for proportional scaling
max_chrom_length = df['query_end'].max()

# Create a plot for each chromosome
fig, axes = plt.subplots(len(chromosomes), 1, figsize=(15, 2.5 * len(chromosomes)), sharex=True)  # Adjusted vertical space

if len(chromosomes) == 1:
    axes = [axes]

# Define a function to calculate density
def calculate_density(chrom_data, chrom_length, window_size):
    density = np.zeros((chrom_length // window_size) + 1)
    for start, end in zip(chrom_data['query_begin'], chrom_data['query_end']):
        start_idx = start // window_size
        end_idx = end // window_size
        for i in range(start_idx, end_idx + 1):
            density[i] += 1
    return density

# Calculate density and max density for each class separately
density_dict = {}
max_density_dict = {cls: 0 for cls in selected_classes}

for chrom in chromosomes:
    chrom_data = df[df['sequence'] == chrom]
    chrom_length = chrom_data['query_end'].max()

    if args.window_size:
        window_size = args.window_size
        density_dict[chrom] = {}
        for cls in selected_classes:
            density = calculate_density(chrom_data[chrom_data['consolidated_class'] == cls], chrom_length, window_size)
            density_dict[chrom][cls] = density
            max_density_dict[cls] = max(max_density_dict[cls], max(density))

for chrom_idx, (ax, chrom) in enumerate(zip(axes, chromosomes)):
    chrom_data = df[df['sequence'] == chrom]
    chrom_length = chrom_data['query_end'].max()

    # Draw the chromosome ideogram as a horizontal bar
    ax.add_patch(mpatches.Rectangle((0, 1.1), chrom_length, 0.1, color='lightgrey', zorder=0))

    if args.window_size:
        window_size = args.window_size
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
            class_data = chrom_data[chrom_data['consolidated_class'] == cls]
            y_position = 1.0 - class_idx * 0.12  # Adjusted spacing between tracks
            for _, row in class_data.iterrows():
                ax.add_patch(mpatches.Rectangle((row['query_begin'], y_position), row['query_end'] - row['query_begin'], 0.1, color=class_colors[row['consolidated_class']], zorder=1))
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
if args.window_size:
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

