#!/usr/bin/env Rscript

# Notes on script
info <- paste("Script to perform a Random Forest model on genetic and phenotypic data.",
              "Adapted from epiGWAS_2025.r on local machine")

# Install and load packages
packages_needed <- c("data.table", "ranger")
packages_missing <- packages_needed[!packages_needed %in% installed.packages()[, "Package"]]
if (length(packages_missing) > 0) install.packages(packages_missing)
invisible(lapply(packages_needed, library, character.only = TRUE))

# Load data
load("CHG_genoPhen.RData")  # loads 'test_2'

# Remove all-NA marker columns (from col 4 onward)
marker_start <- 4
marker_cols <- marker_start:ncol(test_2)
na_cols <- which(colSums(is.na(test_2[, marker_cols])) == nrow(test_2))
if (length(na_cols) > 0) {
  test_2 <- test_2[, -((marker_start - 1) + na_cols)]
  marker_cols <- marker_start:ncol(test_2)  # reassign range
}

# Remove rows where all marker values are NA
na_rows <- which(rowSums(is.na(test_2[, marker_cols])) == length(marker_cols))
if (length(na_rows) > 0) test_2 <- test_2[-na_rows, ]

# Filter rows with missing phenotype or markers
complete_cases <- complete.cases(test_2$blup, test_2[, marker_cols])
rf_data <- data.frame(trait = test_2$blup[complete_cases],
                      test_2[complete_cases, marker_cols])

# Optional: Clean up workspace before modeling
rm(test_2, na_cols, na_rows, complete_cases, marker_cols, marker_start)
gc()

# Fit Random Forest model
rf_model <- ranger(
  trait ~ .,
  data = rf_data,
  num.trees = 500,
  num.threads = parallel::detectCores() - 1,
  importance = "impurity"
)

### Save the model
saveRDS(rf_model, file = "rf_model_CHG.rds") 

# Variable importance table
importance_df <- data.frame(
  marker = names(rf_model$variable.importance),
  importance = rf_model$variable.importance
)
importance_df <- importance_df[order(-importance_df$importance), ]
fwrite(importance_df, "rf_importance_scores.csv")

# Session snapshot
save.image(file = "rf_analysis_workspace.RDat


