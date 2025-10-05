# =============================================================
# Environment File for BLAST Plotting Script
# -------------------------------------------------------------
# This file sets environment variables used by `load_blast_data()`
# in the prep_multiplot_module.R script. It is designed
# to load and preprocess multiple BLAST result files and their
# corresponding chromosome length tables.
#
# FORMAT:
#   - Use `bln_path_N` for each BLAST tabular file (13-column format)
#   - Use `len_N` for the corresponding chromosome length file (2-column TSV)
#   - N should be a matching index (e.g., 1, 2, 3, ...)
#
# USAGE:
#   The script will dynamically detect and load all pairs of
#   `bln_path_N` and `len_N` entries. Each pair represents a
#   dataset to be plotted as a vertical strip of chromosomes.
#
# Note: All file paths should be absolute or relative to the
#       working directory where the script is executed.
# =============================================================

# BLAST file paths
bln_path_1=/path/to/blast_dataset1.bln
bln_path_2=/path/to/blast_dataset2.bln

# Corresponding chromosome length files
len_1=/path/to/chr_lengths1.txt
len_2=/path/to/chr_lengths2.txt
