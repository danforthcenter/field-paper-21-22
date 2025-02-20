# Load  data
load("/path/to/asvTable_noAbnd.Rdata") # Load the asvTable_noAbnd data

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
re <- read_excel("/Users/eflom/Library/CloudStorage/OneDrive-DDPSC/bart_lab/21_22_field_paper/Correlation/spearman_correlation/combined_output_zones.xlsx", sheet = )
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
get_helper_abundance <- function(asv_data, helper_asvs) {
  # Subset only the helper ASVs from the main ASV abundance table
  return(asv_data[, helper_asvs])
}

# Function to find the highest abundance sample for each helper ASV
get_highest_abundance_per_helper <- function(helper_abundance_data) {
  # Apply function to find the sample with highest abundance for each helper ASV
  highest_abundance_per_helper <- apply(helper_abundance_data, 2, function(x) {
    sorted_abundance <- sort(x, decreasing = TRUE)
    top_10_samples <- names(sorted_abundance)[1:10]
    top_10_abundance <- sorted_abundance[1:10]
    c(Sample = top_10_samples, Abundance = top_10_abundance)
  })
  return(data.frame(highest_abundance_per_helper))
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
save_results <- function(highest_abundance_per_helper, highest_combined_sample, output_directory) {
  # Create summary results data frame
  summary_results <- data.frame(
    Sample = highest_abundance_per_helper[1, ],
    Helper_ASV = rep(names(highest_abundance_per_helper[1, ]), each = length(highest_abundance_per_helper[1, ])),
    Abundance = highest_abundance_per_helper[2, ]
  )

  # Save individual helper abundance summary
  write.csv(summary_results, file.path(output_directory, "helper_abundance_summary.csv"), row.names = FALSE)

  # Save highest combined abundance sample info
  write.csv(data.frame(Sample = highest_combined_sample$Sample, Combined_Abundance = highest_combined_sample$Combined_Abundance),
    file.path(output_directory, "combined_helper_abundance.csv"),
    row.names = TRUE
  )

  cat("Results saved to", output_directory, "\n")
}

# Main function to process helper abundance data
process_helper_abundance <- function(asv_data, helper_df, output_directory) {
  # Get the helper ASVs
  helper_asvs <- helper_df$ASV

  # Get helper ASV abundance data from the main table
  helper_abundance_data <- get_helper_abundance(asv_data, helper_asvs)

  # Find the highest abundance sample for each helper ASV
  highest_abundance_per_helper <- get_highest_abundance_per_helper(helper_abundance_data)

  # Get combined abundance across all helper ASVs for each sample
  combined_abundance <- get_combined_helper_abundance(helper_abundance_data)

  # Find the sample with the highest combined abundance
  highest_combined_sample <- get_top_10_combined_abundance_samples(combined_abundance)

  # Save the results to CSV files
  save_results(highest_abundance_per_helper, highest_combined_sample, output_directory)

  cat("Processing complete.\n")
}
