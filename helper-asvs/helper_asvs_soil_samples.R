library(readxl)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)

# Load the asvTable_noAbnd data
asvTable <- read.csv("/Users/eflom/Library/CloudStorage/OneDrive-DDPSC/bart_lab/21_22_field_paper/Abundance/filtered_ASV_noabnd/asvTable_Abnd_25_or_more.csv", row.names = 1)
# Remove rows with zero ASV counts
zero_sum_rows <- rowSums(asvTable) == 0
asvTable <- asvTable[!zero_sum_rows, ]
# Normalize the ASV table to relative abundance
asvTable <- asvTable / rowSums(asvTable) # Convert to relative abundance


# Load the taxa
load("/Users/eflom/Downloads/primerless_taxa_rdp.rdata")
taxa_asvs <- rownames(taxa)

# Load the correlation data for root rhizosphere, root endosphere, and soil
root_rhizosphere <- read_excel("/Users/eflom/Library/CloudStorage/OneDrive-DDPSC/bart_lab/21_22_field_paper/Correlation/spearman_correlation/combined_output_zones.xlsx", sheet = "rr")
root_endosphere <- read_excel("/Users/eflom/Library/CloudStorage/OneDrive-DDPSC/bart_lab/21_22_field_paper/Correlation/spearman_correlation/combined_output_zones.xlsx", sheet = "re")
soil <- read_excel("/Users/eflom/Library/CloudStorage/OneDrive-DDPSC/bart_lab/21_22_field_paper/Correlation/spearman_correlation/combined_output_zones.xlsx", sheet = "s")

# Function to extract the helper ASVs from the main ASV data
get_helper_abundance <- function(asv_data, helper_df) {
  # Arrange by descending Correlation for each data frame using dplyr
  helper_df <- dplyr::arrange(helper_df, desc(Correlation))
  helper_df <- helper_df[1:10, ] # Limit to top 10 helpers
  # Subset only the helper ASVs from the main ASV abundance table
  helper_asvs <- unique(helper_df$ASV)
  valid_helpers <- intersect(helper_asvs, colnames(asv_data))
  # Subset only the valid helper ASVs
  helper_abundance_data <- asv_data[, valid_helpers, drop = FALSE]
  # Check if any valid helpers were found
  if (ncol(helper_abundance_data) == 0) {
    warning("No valid helper ASVs found in the asvTable.")
  }

  return(helper_abundance_data)
}

# Function to find the highest abundance sample for each helper ASV
get_highest_abundance_per_helper <- function(helper_abundance_data) {
  # Initialize empty list to store results
  result_list <- list()

  # Loop through each ASV and extract top 10 samples and abundance
  for (asv in colnames(helper_abundance_data)) {
    # Extract the abundance values for the current ASV
    abundance_values <- helper_abundance_data[[asv]]
    # Get the sorted indices for abundance (in decreasing order)
    sorted_indices <- order(abundance_values, decreasing = TRUE)
    # Sort the abundance values and sample names based on sorted indices
    sorted_abundance <- abundance_values[sorted_indices]
    # Define the number of samples (top 10 or less if not available)
    num_samples <- min(10, length(sorted_abundance)) # Ensure we don't go beyond available samples
    # Sample names (rownames) corresponding to the sorted abundance
    top_10_samples <- rownames(helper_abundance_data)[sorted_indices][1:num_samples]
    top_10_abundance <- sorted_abundance[1:num_samples]
    # Ensure there are at least some non-zero values
    if (length(sorted_abundance) == 0 || all(sorted_abundance == 0)) {
      message(paste("No valid abundance for ASV:", asv))
      return(NULL)
    }

    # Mismatch check
    if (length(top_10_samples) != length(top_10_abundance)) {
      message("Mismatch in top samples and abundance lengths!")
      message("Length of top_10_samples: ", length(top_10_samples))
      message("Length of top_10_abundance: ", length(top_10_abundance))
      return(NULL)
    }

    # Append the results to a list
    result_list[[asv]] <- data.frame(
      Sample = top_10_samples,
      Helper_ASV = rep(asv, num_samples),
      Abundance = top_10_abundance,
      stringsAsFactors = FALSE
    )
  }

  # Combine all ASV results into a single data frame
  result_df <- do.call(rbind, result_list)

  return(result_df)
}


# Function to calculate the combined abundance of all helpers for each sample
get_top_5_samples_for_helper_abundance <- function(helper_abundance_data) {
  # Sum the abundance of all helper ASVs for each sample
  combined_abundance <- rowSums(helper_abundance_data)

  # Get the top 5 samples with the highest combined abundance
  top_5_samples <- names(sort(combined_abundance, decreasing = TRUE)[1:5])
  # Subset the data to include only the top 5 samples
  top_5_samples_data <- helper_abundance_data[top_5_samples, , drop = FALSE]
  return(top_5_samples_data)
}

# Function to save results as CSV
save_results <- function(top_5_samples, helper_df_name, output_directory) {
  # Save highest combined abundance sample info
  write.csv(top_5_samples, file.path(output_directory, paste0(helper_df_name, "_top_5_samples.csv")), row.names = FALSE)
  cat("Results saved to", output_directory, "\n")
}

# Main function to process helper abundance data
process_helper_abundance <- function(asv_data, helper_df, output_directory) {
  # Get helper ASV abundance data from the main table
  helper_abundance_data <- get_helper_abundance(asv_data, helper_df)
  # Find the sample with the highest combined abundance
  top_5_samples <- get_top_5_samples_for_helper_abundance(helper_abundance_data)

  # Reshape the data to long format for ggplot2
  top_5_samples_long <- top_5_samples %>%
    rownames_to_column("Sample") %>%
    pivot_longer(cols = -Sample, names_to = "Helper_ASV", values_to = "Abundance") %>%
    filter(!is.na(Abundance)) # Filter out any rows with NA values (if any)

  # Pull out the helper_df string
  helper_df_name <- as.character(match.call()[[3]])
  # Replace underscores with spaces
  formatted_helper_df_name <- gsub("_", " ", helper_df_name)
  # Capitalize the first letters of each word (using sub to capitalize the first letter of 'root' and 'endosphere')
  formatted_helper_df_name <- sub("(^|\\s)([a-z])", "\\1\\U\\2", formatted_helper_df_name, perl = TRUE)
  formatted_helper_df_name <- sub("(^|\\s)([a-z])", "\\1\\U\\2", formatted_helper_df_name, perl = TRUE)

  # Save the results to CSV files
  save_results(top_5_samples, helper_df_name, output_directory)

  # Save taxa data
  taxa_data <- taxa[taxa_asvs %in% unique(helper_df$ASV), ]
  # Save the merged data to a CSV file
  write.csv(taxa_data, file.path(output_directory, paste0(helper_df_name, "_helper_taxa.csv")), row.names = TRUE)

  ggplot(top_5_samples_long, aes(x = Sample, y = Abundance, fill = Helper_ASV)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(
      title = paste("Top 5 Samples for", formatted_helper_df_name, "Helper ASVs"),
      x = "Samples", y = "Relative Abundance"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(drop = FALSE) # Ensure all Sample levels are shown, even if some have no data

  # Save the plot
  ggsave(
    filename = paste0(output_directory, "/", helper_df_name, "_top_5_samples_abundance_plot.png"),
    width = 12, height = 6
  )

  cat("Processing complete.\n")
}
