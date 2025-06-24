#!/usr/bin/env python3

# --- Setup on commandline ---
awk '/.Gm15/ && $2>=40000000 && $2<=42000000 && $4==91 {print}' $file > sub.glyma.Gm15.91.40-42Mb.tsv
awk '/.Gm15/ && $2>=40000000 && $2<=42000000 && $4==92 {print}' $file > sub.glyma.Gm15.92.40-42Mb.tsv

# --- Setup in Python ---

# Import libraries
import pandas as pd
import numpy as np
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor
from Bio.Align import PairwiseAligner

# Read the file and return the dataframe
def read_file(filename):
    """
    Reads the ULTRA file and returns a dataframe.
    """
    return pd.read_csv(filename, sep="\t", 
                        names=['SequenceName', 'Start', 
                        'Length', 'Period', 'Score', 'Substitutions', 
                        'Insertions', 'Deletions', 'Consensus', 'Sequence'], 
                        usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

df = read_file(filename)
df.description = "Chromosome 15 data from Glycine max"

# 🧬 Your manually chosen representative
representative_seq = 'TGTGAAAAGTTATGACCATTTGAATTTCTCGAGAGCTTCCGTTGTTCAATTTCGAGCGTCTCGATATATTATGCGCCTGAATCGGACATCCG'

# --- Feature extraction ---
def compute_entropy(seq):
    freqs = [seq.count(nuc)/len(seq) for nuc in set(seq)]
    return -sum(f*np.log2(f) for f in freqs if f > 0)

df["GC_Content"] = df["Consensus"].apply(lambda s: (s.count("G") + s.count("C")) / len(s))
df["Entropy"] = df["Consensus"].apply(compute_entropy)
df["Indel_Variability"] = df["Substitutions"] + df["Insertions"] + df["Deletions"]


# Centromere midpoint analysis
centromere_midpoint = 41000000  # Example midpoint for Gm15
df["DistanceFromCentromere"] = abs(df["Start"] - centromere_midpoint)
df["NormalizedDistance"] = df["DistanceFromCentromere"] / df["Length"]


# --- Alignment-based outlier detection ---
aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 1
aligner.mismatch_score = 0
aligner.open_gap_score = 0
aligner.extend_gap_score = 0

# Apply alignment scoring to each row
df["AlignmentScore"] = [
    aligner.score(representative_seq, query_seq) / max(len(representative_seq), len(query_seq))
    for query_seq in df["Consensus"]
]

# Flag sequences that deviate significantly from the representative
df["Align_Outlier"] = df["AlignmentScore"] < 0.90

# --- Feature-based outlier detection ---
features = df[["GC_Content", "Entropy", "Indel_Variability"]].values
df["IsolationForest"] = IsolationForest(contamination=0.2).fit_predict(features)
df["LOF"] = LocalOutlierFactor(n_neighbors=10, contamination=0.2).fit_predict(features)

# Optimizing the outlier detection (particularly for Local outlier factor)
# Check for duplicates in the feature space to optimize neighbors parameter
duplicates = df[["GC_Content", "Entropy", "Indel_Variability"]].duplicated(keep=False)
print(df[duplicates])
# 38 duplicates found in df of 215 rows. If 38 rows are duplicates, a good starting 
# point is n_neighbors ≈ 10–20, or up to ~√N. That is 18% of the data. LOF depends 
# on nearest neighbors, so if many points are exact duplicates, LOF can’t rank their 
# local density meaningfully unless there are enough neighbors to differentiate between them.
# With n_neighbors=2, LOF is effectively comparing each point to 1 or 2 others — not enough when 38+ values are the same.


# Side by side comparison of outlier detection method
# Basic features
features_basic = df[["GC_Content", "Entropy", "Indel_Variability"]].values
df["IF_Basic"] = IsolationForest(contamination=0.2).fit_predict(features_basic)
df["LOF_Basic"] = LocalOutlierFactor(n_neighbors=10, contamination=0.2).fit_predict(features_basic)

# Enhanced features
features_enhanced = df[["GC_Content", "Entropy", "Indel_Variability", "NormalizedDistance", "AlignmentScore"]].values
df["IF_Enhanced"] = IsolationForest(contamination=0.2).fit_predict(features_enhanced)
df["LOF_Enhanced"] = LocalOutlierFactor(n_neighbors=10, contamination=0.2).fit_predict(features_enhanced)

# Count how many outliers were flagged
(df["IF_Basic"] == -1).sum(), (df["IF_Enhanced"] == -1).sum()
# Check how many sequences disagree
(df["IF_Basic"] != df["IF_Enhanced"]).sum()

# --- Visual Comparison (e.g. t-SNE, UMAP, PCA, t-SNE) ---
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt

tsne = TSNE(n_components=2, random_state=42)
embedding1 = tsne.fit_transform(features_basic)
embedding2 = tsne.fit_transform(features_enhanced)

plt.scatter(embedding1[:, 0], embedding1[:, 1], c=(df["IF_Basic"] == -1), cmap="coolwarm")
plt.title("Isolation Forest (Basic Features)")
plt.show()
plt.savefig(IF_basicfeat.png')

plt.clf()

plt.scatter(embedding2[:, 0], embedding2[:, 1], c=(df["IF_Enhanced"] == -1), cmap="coolwarm")
plt.title("Isolation Forest (Enhanced Features)")
plt.show()
plt.savefig(IF_enhancedfeat.png')

# color by other labels too (ie. AlignmentScore)
plt.scatter(embedding1[:, 0], embedding1[:, 1], c=df["AlignmentScore"], cmap="viridis")
# then replicate plotting with IF_Enhanced and LOF_Enhanced

plt.scatter(embedding2[:, 0], embedding2[:, 1], c=df["AlignmentScore"], cmap="viridis")
# then replicate plotting with IF_Enhanced and LOF_Enhanced

# follouw ups
# threshold overlay
low_score = df["AlignmentScore"] < 0.90
plt.scatter(embedding[:, 0], embedding[:, 1], c=low_score, cmap="coolwarm")

# label top matches
top_hits = df["AlignmentScore"] > 0.95
plt.scatter(..., label="Highly similar")

# --- Side-by-Side 2D Scatter Plots of...Isolation Forest vs LOF...before vs after adding features ---
# Dimensionality reduction
embedding = TSNE(n_components=2, random_state=42).fit_transform(...features)

# Create Isolation Forest vs LOF side-by-side plots
fig, axes = plt.subplots(1, 2, figsize=(12, 5))  # 1 row, 2 columns

# Plot 1: Isolation Forest
axes[0].scatter(
    embedding[:, 0], embedding[:, 1],
    c=(df["IsolationForest"] == -1), cmap="coolwarm", s=10
)
axes[0].set_title("Isolation Forest Outliers")

# Plot 2: Local Outlier Factor
axes[1].scatter(
    embedding[:, 0], embedding[:, 1],
    c=(df["LOF"] == -1), cmap="coolwarm", s=10
)
axes[1].set_title("LOF Outliers")

plt.tight_layout()
plt.show()

# Create before vs after adding features side-by-side plots
fig, axes = plt.subplots(2, 2, figsize=(15, 10))

# Plot 1: Basic features Isolation Forest
axes[0,0].scatter(
    embedding1[:, 0], embedding1[:, 1], 
    c=(df["IF_Basic"] == -1), cmap="coolwarm"
)
axes[0,0].set_title("Basic features Isolation Forest")

# Plot 2: Enhanced features Isolation Forest
axes[0,1].scatter(
    embedding2[:, 0], embedding2[:, 1], 
    c=(df["IF_Enhanced"] == -1), cmap="coolwarm"
)
axes[0,1].set_title("Enhanced features Isolation Forest")

# Make side-by-side plots colored by AlignmentScore, DistanceFromCentromere, Forest predictions
# Customize Further with a thrid panel for AlignmentScore as a color gradient

axes[1,0].scatter(
    embedding1[:, 0], embedding1[:, 1],
    c=df["AlignmentScore"], cmap="viridis", s=10
)
axes[1,0].set_title("Basic Alignment Score Gradient")

axes[1,1].scatter(
    embedding2[:, 0], embedding2[:, 1],
    c=df["AlignmentScore"], cmap="viridis", s=10
)
axes[1,1].set_title("Enhanced Alignment Score Gradient")

plt.tight_layout()
plt.show()

# --- View results ---
print(df[[
    "SequenceName", "GC_Content", "Entropy", "Indel_Variability",
    "IsolationForest", "LOF", "AlignmentScore", "Align_Outlier",
    "DistanceFromCentromere", "NormalizedDistance"
]])

# Save the dataframe without outliers and outliers as separate files
# Define your filtering logic
filtered_df = df[
    (df["IF_Enhanced"] == 1) & # choose the feature set you want to filter by
    (df["LOF_Enhanced"] == 1) & # the enhanced features were chosen here
    (df["AlignmentScore"] >= 0.90)
]

# Save to new TSV file
filtered_df.to_csv("filtered_repeats.tsv", sep="\t", index=False)

print(f"{len(filtered_df)} rows retained out of {len(df)}")

outliers = df[~df.index.isin(filtered_df.index)] 
outliers.to_csv("removed_outliers.tsv", sep="\t", index=False)

# breakdown report
print("Total rows:", len(df))
print("Filtered rows retained:", len(filtered_df))
print("Removed as IsolationForest outlier:", (df["IsolationForest"] == -1).sum())
print("Removed as LOF outlier:", (df["LOF"] == -1).sum())
print("Removed for low alignment (< 0.90):", (df["AlignmentScore"] < 0.90).sum())


# --- On the command line now...create FASTA files to check the alignments visually after filtering ---
for filename in sub*tsv; do
    ref=$(echo $filename | perl -pe 's/^sub\.(.+)\.tsv/$1/')
    echo "$ref"
    awk '{print ">"$1"_"$2"_"$3"_"$4"\n"$9}' $filename > sub.$ref.fna
done   
