#' load_blast_data
#'
#' Loads and prepares multiple BLAST datasets and their associated chromosome lengths
#' for comparative plotting. File paths are defined in a user-supplied `.env_r` file.
#'
#' The `.env_r` file must define variables named in the form:
#'   - `bln_path_1`, `bln_path_2`, ..., for each BLAST file to be loaded
#'   - `len_1`, `len_2`, ..., for each corresponding chromosome length file
#'
#' The number of datasets is inferred by counting how many `bln_path_N` entries exist.
#'
#' @param env_file Path to a `.env_r` file containing named paths to input data
#'
#' @return A list containing:
#'   - blast_list: a list of cleaned, merged BLAST data.frames (one per dataset)
#'   - max_chr_length: maximum chromosome length across all datasets
#'
#' @examples
#'   result <- load_blast_data(".env_r")
#'   result$blast_list[[1]]  # first dataset
load_blast_data <- function(env_file) {
  readRenviron(env_file)

  # Gather all environment variables starting with "bln_path_"
  env_vars <- Sys.getenv()
  bln_paths <- env_vars[grep("^bln_path_\\d+$", names(env_vars))]
  len_paths <- env_vars[grep("^len_\\d+$", names(env_vars))]

  if (length(bln_paths) != length(len_paths)) {
    stop("Mismatch: each 'bln_path_N' must have a matching 'len_N' entry.")
  }

  blast_list <- list()
  all_chr_lengths <- c()

  for (i in seq_along(bln_paths)) {
    # Read blast file
    blast <- read.table(bln_paths[i], header = FALSE, sep = "\t", col.names = c(
      "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
      "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sseq"
    ))

    # Read chr lengths
    chr_len <- read.table(len_paths[i], header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(chr_len) <- c("sseqid", "chr_length")

    # Simplify chromosome names
    blast$sseqid <- sub(".*\\.", "", blast$sseqid)
    chr_len$sseqid <- sub(".*\\.", "", chr_len$sseqid)

    # Merge and sort
    blast <- merge(blast, chr_len, by = "sseqid")
    chrom_nums <- as.numeric(gsub("Chr", "", sort(unique(chr_len$sseqid))))
    chrom_order <- paste0("Chr", sprintf("%02d", sort(chrom_nums)))
    blast$sseqid <- factor(blast$sseqid, levels = chrom_order)

    blast_list[[i]] <- blast
    all_chr_lengths <- c(all_chr_lengths, blast$chr_length)
  }

  max_chr_length <- max(all_chr_lengths, na.rm = TRUE)

  return(list(
    blast_list = blast_list,
    max_chr_length = max_chr_length
  ))
}


# Define plotting colors
#custom_colors <- c("#C8102E", "#F1BE48") # ISU colors (Cardinal and Gold)
# other color options:
#custom_colors <- c("#1B9E77", "#4C5C68", "#D95F02")  # Clean Genomics set
# or
#custom_colors <- c("#B22222", "#4B0082", "#E69F00")  # Bold Academic
# or
custom_colors <- c("#556B2F", "#8B4000", "#6C7A89")  # Earthy Modern

