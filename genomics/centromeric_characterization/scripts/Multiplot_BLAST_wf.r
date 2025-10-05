# Multiplot of BLAST data

# Set your working directory
setwd("C:/Users/bdjor/Desktop/temp/genome/blastout/") 

#source modules for load_blast_data(), plot_blast(), and plot_blast_grid()
source("scripts/prep_multiplot_module.R")
source("scripts/multiplot__blast_module.R")

# parse data for passing to plotting function
bl  <- prep$blast_list
ids <- names(bl)
mx  <- prep$max_chr_length

# set custom titles
alt <- c("Soybean v6", "Wild Soybean", "Species C", "Species D")

# create the plot grid
final_plot <- plot_blast_grid(
    blast_list     = prep$blast_list,
    ids            = names(prep$blast_list),
    max_chr_length = prep$max_chr_length,
    colors         = custom_colors, # or use default colors by leaving NULL
    alt_titles     = alt,           # or use NULL to just use ids
    title          = "Locations of Putative X Centromeric Repeat",
    fill_method = "wheel" # or "custom"
)

# save the plot
ggsave("patchwork_multi_blasts.pdf", plot = final_plot, width = 14, height = 14)


# Define plotting colors
#custom_colors <- c("#C8102E", "#F1BE48") # ISU colors (Cardinal and Gold)
# other color options:
#custom_colors <- c("#1B9E77", "#4C5C68", "#D95F02")  # Clean Genomics set
# or
#custom_colors <- c("#B22222", "#4B0082", "#E69F00")  # Bold Academic
# or
custom_colors <- c("#556B2F", "#8B4000", "#6C7A89")  # Earthy Modern

