#!/usr/bin/env Rscript

###### Author: Brandon D. Jordan
###### Started 07/08/2024
###### Obj: Visualization tools from the R for centromeric repeat research.
######      The scripts provided here were useful for data curation and analysis of candidate centromeric repeats
######      Commandline scripts are provided where suffcient to assist data preparation for plotting.

# Create ridgeline multi-plots using ggplot2. This may represent an advanced version of the traditional histogram.

    # Data prep for plotting all repeats in a multiplot of all or a subset of chromosomes
    
    for i in {01..11}; do echo Chr$i ; awk ''/'Chr'$i'/ {print int($2/1000000)}' $file | sort -n | uniq -c | awk '{print $0" "'$i'}' | sed 's/^ *//' >> test.1.txt; done

    # Data prep for plotting a select repeat or range in a multiplot of all or a subset of chromosomes (note this does not prepare for resolving between repeat lengths in the plot)
    
    for i in {01..11}; do echo Chr$i ; awk ''/'Chr'$i'/ && $4==193 {print int($2/1000000)}' $file | sort -n | uniq -c | awk '{print $0" "'$i'}' | sed 's/^ *//' >> test.2.txt; done

    # Data prep for plotting select repeats in a multiplot of all or a subset of chromosomes
    
    i=1 
    for x in 71 144; do echo Chr$i $x; awk ''/'Chr'$i'/ && $4=='$x' {print int($2/1000000)}' $file | sort -n | uniq -c | awk '{print $0" "'$x'}' | sed 's/^ *//' >> test2.txt; done

    # Data prep for plotting select repeats in individual plots for all or a subset of chromosomes

    for i in {01..11}; do for x in 84 193; do echo Chr$i $x; awk ''/'Chr'$i'/ && $4=='$x' {print int($2/1000000)}' $file | sort -n | uniq -c | awk '{print $0" "'$x'}' | sed 's/^ *//' >> test.chr$i.txt; done; done

Based on common background repeat types in plant genomes, and avoiding known centromeric sizes (170–200 bp range), I selecting the most
common tandem repeats of ~30–60 bp to serve as background for centromeric signals. Shorter satellite repeats (~30–60 bp) often occur in subtelomeric and intergenic regions and are less likely to be enriched at centromeres. These sizes are frequently found in RepeatExplorer outputs and tend to be more evenly
distributed in plants, including legumes.

### Import libraries

    library(ggplot2)
    library(ggridges)
    library(viridis)
    library(dplyr)

#### Move to working directory

    getwd()
    setwd("C:/.../.../")

#### Import files

    chrom_data <- read.table('test.txt', header = FALSE, sep = " ", col.names = c("frequency", "position", "chromosome")) # the 3rd column will change according to the data. It is used for distinguishing between plots.
    chrom_data$chromosome <- factor(chrom_data$chromosome)

#### Plot data for all or a select repeat across multpile chromosomes with geom_density_ridges_gradient()

    plot <- ggplot(chrom_data, aes(x = position, y = chromosome, height = frequency, fill = after_stat(x), group = chromosome)) +
    geom_density_ridges_gradient(stat = "identity", scale = 2, rel_min_height = 0.01) +
    scale_fill_viridis_c(name = "Frequency", option = "C") +
    theme_minimal() +
    labs(title = "Apipr Whole Genome Distribution of Tandem Repeats",
        x = "Genomic Position",
        y = "Chromosome")
    print(plot)
    ggsave("Apipr_chromosome_ridgeline_plot.png", plot = plot, width = 10, height = 6, dpi = 300)

#### Plot data for all or a select repeat across multpile chromosomes with geom_density_ridges_gradient()

    plot <- ggplot(chrom_data, aes(x = position, y = chromosome, height = frequency, fill = after_stat(x), group = chromosome)) +
    geom_density_ridges_gradient(stat = "identity", scale = 2, rel_min_height = 0.01) +
    scale_fill_viridis_c(name = "Frequency", option = "C") +
    theme_minimal() +
    labs(title = "Apiam Whole Genome Distribution of 193 bp Tandem Repeats",
        x = "Genomic Position",
        y = "Chromosome")
    print(plot)
    ggsave("Apiam_193_chromosome_ridgeline_plot.png", plot = plot, width = 10, height = 6, dpi = 300)

#### Plot data for a single chromosome with geom_density_ridges_gradient()

    chrom_data$repeat. <- factor(chrom_data$repeat.)
    plot <- ggplot(chrom_data, aes(x = position, y = repeat., height = frequency, fill = after_stat(x), group = repeat.)) +
    geom_density_ridges_gradient(stat = "identity", scale = 2, rel_min_height = 0.01) +
    scale_fill_viridis_c(name = "Frequency", option = "C") +
    theme_minimal() +
    labs(title = "Chr.1 Frequency Ridgeline Plot (Repeat size: 104, 71)",
        x = "Genomic Position",
        y = "Repeat")
    print(plot)
    ggsave("Chr1_ridgeline_plot.png", plot = plot, width = 10, height = 6, dpi = 300)


#### My favorite for plotting several repeats on individual chromosomes geom_ridgeline()

    chrom_data <- read.table('test1.txt', header = FALSE, sep = " ", col.names = c("frequency", "position", "repeat"))

##### Convert the repeat column to a factor

    repeat_order <- chrom_data %>%
    group_by(repeat.) %>%
    summarise(max_freq = max(frequency)) %>%
    arrange(max_freq) %>%   # highest freq first (goes to back)
    pull(repeat.)

    chrom_data$repeat. <- factor(chrom_data$repeat., levels = repeat_order)
    plot <- ggplot(chrom_data, aes(x = position, y = repeat., height = frequency, fill = repeat., group = repeat.)) +
    geom_ridgeline(color = "black", alpha = 0.8) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal() +
    labs(fill = "Repeat Group")
    ggsave("Cerchi_chr1_ridgeline_plot.png", plot = plot, width = 10, height = 6, dpi = 300)

#### Same methodolgy basically for plotting batches of individual chr plots, edit input and output filenames accordingly

    for (x in 1:11) {
    # Create input and output filenames
    chr_str <- sprintf("%02d", x)  # ensures two-digit format
    input_file <- paste0("test.2.chr", chr_str, ".36x84x193.txt") # EDIT HERE
    output_file <- paste0("Apipr_chr", chr_str, "_ridgeline_plot.png") # EDIT HERE

    # Read the data
    chrom_data <- read.table(input_file, header = FALSE, sep = " ",
                            col.names = c("frequency", "position", "repeat"))

    # Convert the repeat column to a factor
    repeat_order <- chrom_data %>%
    group_by(repeat.) %>%
    summarise(max_freq = max(frequency)) %>%
    arrange(max_freq) %>%    # highest freq first (goes to back)
    pull(repeat.)

    chrom_data$repeat. <- factor(chrom_data$repeat., levels = repeat_order)

    # Build the plot
    plot <- ggplot(chrom_data, aes(x = position, y = repeat., height = frequency,
                                    fill = repeat., group = repeat.)) +
        geom_ridgeline(color = "black", alpha = 0.8) +
        scale_fill_brewer(palette = "Set2") +
        theme_minimal() +
        labs(fill = "Repeat Group",
            title = paste("Ridgeline Plot - Chr", chr_str))

    # Save the plot
    ggsave(output_file, plot = plot, width = 10, height = 6, dpi = 300)

    cat("Saved plot:", output_file, "\n")
    }

# If you really want a gradient fill based on position, then use a continuous scale like scale_fill_viridis_c():

    plot <- ggplot(chrom_data, aes(x = position, y = repeat., height = frequency, fill = after_stat(x), group = repeat.)) +
    geom_ridgeline(color = "black", alpha = 0.8) +
    scale_fill_viridis_c(option = "C", name = "Position") +
    theme_minimal()

    print(plot)

## Other Palette Options
## Try any of these in scale_fill_brewer(palette = ...):

    "Set1", "Set2", "Set3" — qualitative colors

    "Blues", "Greens" — for gradients (if your y is numeric)

## Or use scale_fill_manual() to set custom colors:

    scale_fill_manual(values = c("red", "blue", "green"))
