#!/usr/bin/env Rscript
# Load necessary libraries
library(readr)
library(ggplot2)
library(tidyverse)
library(Hmisc)
library(lintr)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Error: Incorrect number of arguments. Expected 3: BLO_consensus_file, asv_data_file, output_directory.")
}

blo_consensus_file <- args[1]
asv_data_file <- args[2]
output_directory <- args[3]
print(args)
# Read BLO_consensus_file and extract the second column (sseqid)
asv_list <- tryCatch(
  {
    read_tsv(BLO_consensus_file, col_names = TRUE, show_col_types = FALSE) %>%
      pull(2) %>% # Extract the second column (sseqid)
      unique() %>% # Ensure unique ASVs
      na.omit() # Remove NAs
  },
  error = function(e) {
    stop("Error reading BLO_consensus_file: ", e)
  }
)

# Ensure the ASV list is valid
if (length(asv_list) == 0) {
  stop("Error: No valid ASVs found in BLO_consensus_file.")
}

# Read ASV abundance table
asv_data_numeric <- tryCatch(
  {
    as.data.frame(read_csv(asv_data_file, show_col_types = FALSE))
  },
  error = function(e) {
    stop("Error reading asv_data_file: ", e)
  }
)

rownames(asv_data_numeric) <- asv_data_numeric[[1]]
asv_data_numeric <- asv_data_numeric[, -1]

# Function to analyze ASV correlation
analyze_asv_correlation <- function(ASV_of_interest, base_save_path, top_n = 10) {
  # Check if ASV exists in the dataset
  if (!(ASV_of_interest %in% colnames(asv_data_numeric))) {
    message(paste("ASV", ASV_of_interest, "not found in the dataset"))
    return(NULL)
  }

  # Identify samples where ASV_of_interest is non-zero
  asv_of_interest_samples <- asv_data_numeric[, ASV_of_interest] != 0

  overlapping_asvs <- colnames(asv_data_numeric)[colSums(asv_of_interest_samples & (asv_data_numeric != 0)) > 0]
  overlapping_asvs <- setdiff(overlapping_asvs, ASV_of_interest) # Remove the ASV itself

  # Skip ASV if no overlapping data is found
  if (length(overlapping_asvs) == 0) {
    message(paste("Skipping ASV", ASV_of_interest, "- No overlapping data found"))
    return(NULL)
  }

  # Define save path
  save_path <- file.path(base_save_path, ASV_of_interest)
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }

  # Subset data to only include overlapping ASVs
  asv_data_with_overlaps <- asv_data_numeric[, overlapping_asvs, drop = FALSE]


  # Initialize vectors for correlation and p-values
  cor_values <- numeric(length(overlapping_asvs))
  p_values <- numeric(length(overlapping_asvs))
  names(cor_values) <- overlapping_asvs
  names(p_values) <- overlapping_asvs

  # Perform Spearman correlation
  for (current_asv in overlapping_asvs) {
    if (var(asv_data_with_overlaps[, current_asv]) > 0) { # Ensure variation
      correlation_result <- cor.test(asv_data_numeric[, ASV_of_interest],
        asv_data_with_overlaps[, current_asv],
        method = "spearman", exact = FALSE
      )
      cor_values[current_asv] <- correlation_result$estimate
      p_values[current_asv] <- correlation_result$p.value
    } else {
      cor_values[current_asv] <- NA
      p_values[current_asv] <- NA
    }
  }

  # Create results dataframe and filter by p-value < 0.05
  cor_results <- data.frame(
    ASV = overlapping_asvs,
    Correlation = cor_values,
    P_Value = p_values
  ) %>%
    arrange(P_Value) %>% # Sort by lowest p-value
    filter(P_Value < 0.05, !is.na(Correlation))

  # Save significant results
  if (nrow(cor_results) > 0) {
    write_csv(cor_results, file.path(save_path, paste0("correlation_results_", ASV_of_interest, ".csv")))

    # Generate abundance plots for top N correlated ASVs
    top_results <- head(cor_results, top_n)

    for (asv_correlated in top_results$ASV) {
      plot_file <- file.path(save_path, paste0("abundance_plot_", asv_correlated, ".png"))

      p <- ggplot(asv_data_numeric, aes(x = .data[[ASV_of_interest]], y = .data[[asv_correlated]])) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "lm", se = FALSE, color = "blue") +
        labs(
          title = paste("Abundance correlation:", ASV_of_interest, "vs", asv_correlated),
          x = ASV_of_interest, y = asv_correlated
        ) +
        theme_minimal()

      ggsave(plot_file, plot = p, width = 6, height = 4, dpi = 300)
    }
  }
}

# Loop through ASVs and analyze them
for (ASV in asv_list) {
  analyze_asv_correlation(ASV, output_directory)
}

message("Processing complete.")
