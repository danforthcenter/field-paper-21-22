# Load  data
load("/Users/eflom/Library/CloudStorage/OneDrive-DDPSC/bart_lab/21_22_field_paper/asvTable_noAbnd.rdata")
# Load the asvTable_noAbnd data
asvTable <- read.csv("/Users/eflom/Library/CloudStorage/OneDrive-DDPSC/bart_lab/21_22_field_paper/Abundance/filtered_ASV_noabnd/asvTable_Abnd_25_or_more.csv", row.names = 1)
asvTable <- asvTable / rowSums(asvTable) # Convert to relative abundance

# Load the taxa
load("/Users/eflom/Downloads/primerless_taxa_rdp.rdata")
taxa_asvs <- rownames(taxa)

library("readxl")
sheet_names <- excel_sheets(
  "/Users/eflom/Library/CloudStorage/OneDrive-DDPSC/bart_lab/21_22_field_paper/Correlation/spearman_correlation/combined_output_zones.xlsx"
)
df_list <- lapply(sheet_names, function(sheet) {
  read_excel(
    "/Users/eflom/Library/CloudStorage/OneDrive-DDPSC/bart_lab/21_22_field_paper/Correlation/spearman_correlation/combined_output_zones.xlsx",
    sheet = sheet
  )
})
names(df_list) <- sheet_names
rr <- read_excel("/Users/eflom/Library/CloudStorage/OneDrive-DDPSC/bart_lab/21_22_field_paper/Correlation/spearman_correlation/combined_output_zones.xlsx", sheet = "rr")
re <- read_excel("/Users/eflom/Library/CloudStorage/OneDrive-DDPSC/bart_lab/21_22_field_paper/Correlation/spearman_correlation/combined_output_zones.xlsx", sheet = "re")
s <- read_excel("/Users/eflom/Library/CloudStorage/OneDrive-DDPSC/bart_lab/21_22_field_paper/Correlation/spearman_correlation/combined_output_zones.xlsx", sheet = "s")

# Extract the "ASV" column (helper ASVs) from the tibble
s_helper_asvs <- unique(s$ASV)
rr_helper_asvs <- unique(rr$ASV)
re_helper_asvs <- unique(re$ASV)

# Function to pull out samples with at least one helper ASV
get_samples_with_helpers <- function(helper_asvs, asv_data_numeric) {
  # Check if helper ASVs exist in the dataset
  missing_asvs <- setdiff(helper_asvs, colnames(asv_data_numeric))
  if (length(missing_asvs) > 0) {
    warning(paste("Missing helper ASVs in the dataset:", paste(missing_asvs, collapse = ", ")))
  }

  # Subset the data where at least one of the helper ASVs has a non-zero value
  helper_columns <- colnames(asv_data_numeric)[colnames(asv_data_numeric) %in% helper_asvs]
  samples_with_helpers <- asv_data_numeric[rowSums(asv_data_numeric[, helper_columns] > 0) > 0, ]

  # Return the samples with helper ASVs
  return(samples_with_helpers)
}

library(dplyr)
# Function to extract the helper ASVs from the main ASV data
get_helper_abundance <- function(asv_data, helper_df) {
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
get_combined_helper_abundance <- function(helper_abundance_data) {
  # Sum the abundance of all helper ASVs for each sample
  combined_abundance <- rowSums(helper_abundance_data)
  return(combined_abundance)
}

# Function to get the top 10 samples with highest combined abundance
get_top_10_combined_abundance_samples <- function(combined_abundance) {
  sorted_combined <- sort(combined_abundance, decreasing = TRUE)
  top_10_samples <- names(sorted_combined)[1:10]
  top_10_combined_abundance <- sorted_combined[1:10]
  list(Sample = top_10_samples, Combined_Abundance = top_10_combined_abundance)
}

# Function to save results as CSV
save_results <- function(highest_abundance_per_helper, highest_combined_sample, helper_df_name, output_directory) {
  # Ensure highest_abundance_per_helper is formatted correctly
  summary_results <- data.frame(
    Sample = as.character(highest_abundance_per_helper$Sample), # Ensure it's character
    Helper_ASV = highest_abundance_per_helper$Helper_ASV,
    Abundance = as.numeric(highest_abundance_per_helper$Abundance) # Ensure numeric
  )

  # Save individual helper abundance summary
  write.csv(summary_results, file.path(output_directory, paste0(helper_df_name, "_helper_abundance_summary.csv")), row.names = FALSE)

  # Save highest combined abundance sample info
  write.csv(
    data.frame(
      Sample = highest_combined_sample$Sample,
      Combined_Abundance = highest_combined_sample$Combined_Abundance
    ),
    file.path(output_directory, paste0(helper_df_name, "_combined_helper_abundance.csv")),
    row.names = FALSE
  )

  cat("Results saved to", output_directory, "\n")
}

# Main function to process helper abundance data
process_helper_abundance <- function(asv_data, helper_df, output_directory) {
  # Get helper ASV abundance data from the main table
  helper_abundance_data <- get_helper_abundance(asv_data, helper_df)
  # Find the highest abundance sample for each helper ASV
  highest_abundance_per_helper <- get_highest_abundance_per_helper(helper_abundance_data)
  # Get combined abundance across all helper ASVs for each sample
  combined_abundance <- get_combined_helper_abundance(helper_abundance_data)
  # Find the sample with the highest combined abundance
  highest_combined_sample <- get_top_10_combined_abundance_samples(combined_abundance)
  # Pull out helper_df character name
  helper_df_name <- as.character(match.call()[[3]])
  # Save the results to CSV files
  save_results(highest_abundance_per_helper, highest_combined_sample, helper_df_name, output_directory)

  taxa_data <- taxa[taxa_asvs %in% unique(helper_df$ASV), ]


  # Save the merged data to a CSV file
  write.csv(taxa_data, file.path(output_directory, paste0(helper_df_name, "_helper_taxa.csv")), row.names = TRUE)

  cat("Processing complete.\n")
}
