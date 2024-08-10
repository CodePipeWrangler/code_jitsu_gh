import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse
import matplotlib.colors as mcolors

def parse_gff(file_path):
    with open(file_path, 'r') as file:
        lines = [line.strip() for line in file if not line.startswith('##')]
    data = [line.split('\t') for line in lines if 'scaffold' not in line.split('\t')[0]]
    columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    df = pd.DataFrame(data, columns=columns)
    df = df[df['type'] == 'dispersed_repeat']
    df['start'] = pd.to_numeric(df['start'])
    df['end'] = pd.to_numeric(df['end'])
    return df

def count_genes_in_ranges(df, window_size):
    gene_counts = {}
    df['sequence_id'] = df['seqid']

    for _, row in df.iterrows():
        start_position = row['start']
        end_position = row['end']
        sequence_id = row['sequence_id']

        start_range = (start_position - 1) // window_size * window_size + 1
        end_range = start_range + window_size - 1

        gene_counts.setdefault(sequence_id, {})
        gene_counts[sequence_id].setdefault((start_range, end_range), 0)
        gene_counts[sequence_id][(start_range, end_range)] += 1

    return gene_counts

def get_colormap(name):
    if name == 'red':
        return mcolors.LinearSegmentedColormap.from_list("", ["white", "red"])
    elif name == 'green':
        return mcolors.LinearSegmentedColormap.from_list("", ["white", "green"])
    elif name == 'blue':
        return mcolors.LinearSegmentedColormap.from_list("", ["white", "blue"])
    elif name == 'yellow':
        return mcolors.LinearSegmentedColormap.from_list("", ["white", "yellow"])
    elif name == 'red-green':
        return mcolors.LinearSegmentedColormap.from_list("", ["red", "green"])
    elif name == 'blue-yellow':
        return mcolors.LinearSegmentedColormap.from_list("", ["blue", "yellow"])
    else:
        return plt.cm.viridis

# Set up argument parser
parser = argparse.ArgumentParser(description='Generate a plot with gene counts in specified ranges.')
parser.add_argument('-i', '--input', type=str, required=True, help='Path to the input file')
parser.add_argument('-w', '--window', type=int, required=True, help='Window size for counting genes')
parser.add_argument('-c', '--color', type=str, choices=['red', 'green', 'blue', 'yellow', 'red-green', 'blue-yellow'], default='viridis', help='Color scheme for the plot')
args = parser.parse_args()

# Read and process the input file
df = parse_gff(args.input)

# Count genes in specified ranges
gene_counts = count_genes_in_ranges(df, args.window)

# Convert gene counts to DataFrame
rows = []
for sequence_id, ranges in gene_counts.items():
    max_gene_count = max(ranges.values())
    for (start_range, end_range), gene_count in sorted(ranges.items()):
        percentage = (gene_count / max_gene_count) if max_gene_count > 0 else 0
        rows.append([sequence_id, start_range, end_range, gene_count, percentage])

result_df = pd.DataFrame(rows, columns=["Sequence ID", "Start Position", "End Position", "Number of Genes", "Percentage"])

# Plotting
sequences = result_df['Sequence ID'].unique()
max_chrom_end = df['end'].max()  # Find the maximum chromosome length for scaling

fig, axes = plt.subplots(len(sequences), 1, figsize=(15, 2 * len(sequences)), sharex=True, constrained_layout=True)

if len(sequences) == 1:
    axes = [axes]

colormap = get_colormap(args.color)

# Debugging: Print the maximum end for each chromosome
chrom_lengths = df.groupby('seqid')['end'].max().to_dict()
#print("Chromosome lengths:", chrom_lengths)

for ax, seq in zip(axes, sequences):
    seq_data = result_df[result_df['Sequence ID'] == seq]
    max_gene_count = seq_data['Number of Genes'].max()
    chrom_end = chrom_lengths[seq]  # Use the correct maximum end for the current chromosome

    ax.add_patch(mpatches.Rectangle((0, 0.4), chrom_end, 0.2, color='lightgrey', zorder=0))

    for _, row in seq_data.iterrows():
        color = colormap(row['Number of Genes'] / max_gene_count)
        ax.add_patch(mpatches.Rectangle((row['Start Position'], 0.4), row['End Position'] - row['Start Position'], 0.2, color=color, zorder=1))

    ax.set_xlim(0, max_chrom_end)  # Set the same x-axis limit for all subplots
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_title(f'Sequence ID: {seq}', fontsize=10)

norm = plt.Normalize(0, max(result_df['Number of Genes']))
sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
sm.set_array([])

# Add colorbar
fig.colorbar(sm, ax=axes, orientation='vertical', label='Number of Repeats', shrink=0.25)

output_file = 'Repeat_density_plot.svg'
try:
    plt.savefig(output_file)
    print(f"Plot saved to {output_file}")
except Exception as e:
    print(f"Error saving the file: {e}")
    raise

plt.show()

