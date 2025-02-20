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

  valid_helpers <- intersect(helper_asvs$ASV, colnames(asv_data))
  # Subset only the valid helper ASVs
  helper_abundance_data <- asv_data[, valid_helpers, drop = FALSE]
  return(helper_abundance_data)
}

# Function to find the highest abundance sample for each helper ASV
get_highest_abundance_per_helper <- function(helper_abundance_data) {
  # Initialize empty vectors to store the results
  sample_names <- c()
  abundance_values <- c()
  helper_asvs <- c()

  # Loop through each ASV and extract top 10 samples and abundance
  lapply(colnames(helper_abundance_data), function(asv) {
    # Initialize the vectors in the global environment
    sample_names <<- c()
    abundance_values <<- c()
    helper_asvs <<- c()
    # Extract the abundance values for the current ASV
    abundance_values <- helper_abundance_data[[asv]]

    # Get the sorted indices for abundance (in decreasing order)
    sorted_indices <- order(abundance_values, decreasing = TRUE)

    # Sort the abundance values and sample names based on sorted indices
    sorted_abundance <- abundance_values[sorted_indices]

    # Sample names (rownames) corresponding to the sorted abundance
    top_10_samples <- rownames(helper_abundance_data)[sorted_indices][1:10]
    top_10_abundance <- sorted_abundance[1:10] # Sorted abundance values

    # Debugging output
    message(paste("Checking ASV:", asv))
    message("Sorted abundance: ", paste(sorted_abundance, collapse = ", "))

    # Ensure there are at least some non-zero values
    if (length(sorted_abundance) == 0 || all(sorted_abundance == 0)) {
      message(paste("No valid abundance for ASV:", asv))
      return(NULL)
    }

    # Mismatch check
    if (length(top_10_samples) != length(top_10_abundance)) {
      message("Mismatch in top 10 sample and abundance lengths!")
      message("Length of top_10_samples: ", length(top_10_samples))
      message("Length of top_10_abundance: ", length(top_10_abundance))
      return(NULL)
    }

    # Append the results directly to the vectors
    sample_names <<- c(sample_names, top_10_samples)
    abundance_values <<- c(abundance_values, top_10_abundance)
    helper_asvs <<- c(helper_asvs, rep(asv, length(top_10_samples)))
  })

  # Combine the vectors into a data frame
  result_df <- data.frame(
    Sample = sample_names,
    Helper_ASV = helper_asvs,
    Abundance = abundance_values,
    stringsAsFactors = FALSE
  )

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
save_results <- function(highest_abundance_per_helper, highest_combined_sample, output_directory) {
  # Ensure highest_abundance_per_helper is formatted correctly
  summary_results <- data.frame(
    Sample = as.character(highest_abundance_per_helper[1, ]), # Ensure it's character
    Helper_ASV = rep(colnames(highest_abundance_per_helper), each = nrow(highest_abundance_per_helper) / 2),
    Abundance = as.numeric(highest_abundance_per_helper[2, ]) # Ensure numeric
  )

  # Save individual helper abundance summary
  write.csv(summary_results, file.path(output_directory, "helper_abundance_summary.csv"), row.names = FALSE)

  # Save highest combined abundance sample info
  write.csv(
    data.frame(
      Sample = highest_combined_sample$Sample,
      Combined_Abundance = highest_combined_sample$Combined_Abundance
    ),
    file.path(output_directory, "combined_helper_abundance.csv"),
    row.names = FALSE
  )

  cat("Results saved to", output_directory, "\n")
}

# Main function to process helper abundance data
process_helper_abundance <- function(asv_data, helper_df, output_directory) {
  # Get the helper ASVs
  helper_asvs <- unique(helper_df$ASV)

  # Get helper ASV abundance data from the main table
  helper_abundance_data <- get_helper_abundance(asv_data, helper_df)
  print(helper_abundance_data)
  # Find the highest abundance sample for each helper ASV
  highest_abundance_per_helper <- get_highest_abundance_per_helper(helper_abundance_data)
  print(highest_abundance_per_helper)
  # Get combined abundance across all helper ASVs for each sample
  combined_abundance <- get_combined_helper_abundance(helper_abundance_data)
  print(combined_abundance)
  # Find the sample with the highest combined abundance
  highest_combined_sample <- get_top_10_combined_abundance_samples(combined_abundance)
  print(highest_combined_sample)
  # Save the results to CSV files
  save_results(highest_abundance_per_helper, highest_combined_sample, output_directory)

  cat("Processing complete.\n")
}
