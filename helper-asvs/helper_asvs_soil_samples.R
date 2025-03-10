library(readxl)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

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
  n_helpers <- min(10, nrow(helper_df)) # Get the smaller value between 10 and the number of helpers
  helper_df <- helper_df[1:n_helpers, ] # Take top N (10) helpers
  # Subset only the helper ASVs from the main ASV abundance table
  helper_asvs <- unique(helper_df$ASV)
  valid_helpers <- intersect(helper_asvs, colnames(asv_data))
  # Subset only the valid helper ASVs
  helper_abundance_data <- asv_data[rowSums(asv_data[, valid_helpers, drop = FALSE]) > 0, valid_helpers, drop = FALSE]
  # Check if any valid helpers were found
  if (ncol(helper_abundance_data) == 0) {
    warning("No valid helper ASVs found in the asvTable.")
  }

  return(helper_abundance_data)
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
  # Get the top 5 samples with the highest combined abundance
  top_5_samples <- get_top_5_samples_for_helper_abundance(helper_abundance_data)

  ## ---- Top 5 Samples ----
  # Reshape the data to long format for ggplot2
  top_5_samples_long <- top_5_samples %>%
    rownames_to_column("Sample") %>%
    pivot_longer(cols = -Sample, names_to = "Helper_ASV", values_to = "Abundance") %>%
    filter(!is.na(Abundance)) # Filter out any rows with NA values (if any)
  # Add a new column for the helper ASV label with the correlation coefficient
  top_5_samples_long <- top_5_samples_long %>%
    mutate(
      Correlation = helper_df$Correlation[match(Helper_ASV, helper_df$ASV)], # Match correlation values
      Helper_Label = paste(Helper_ASV, "(", round(Correlation, 2), ")", sep = "")
    )
  ## ---- All samples ----
  # Convert the full helper abundance data to long format for ggplot2
  all_samples_long <- helper_abundance_data %>%
    rownames_to_column("Sample") %>%
    pivot_longer(cols = -Sample, names_to = "Helper_ASV", values_to = "Abundance") %>%
    filter(!is.na(Abundance)) # Remove NA values

  # Add a new column for the helper ASV label with the correlation coefficient
  all_samples_long <- all_samples_long %>%
    mutate(
      Helper_Label = paste(Helper_ASV, " (",
        round(helper_df$Correlation[match(Helper_ASV, helper_df$ASV)], 2),
        ")",
        sep = ""
      )
    )
  # Pull out the helper_df string
  helper_df_name <- as.character(match.call()[[3]])
  formatted_helper_df_name <- gsub("_", " ", helper_df_name) # Replace underscores with spaces
  # Capitalize the first letters of each word (using sub to capitalize the first letter of 'root' and 'endosphere')
  formatted_helper_df_name <- sub("(^|\\s)([a-z])", "\\1\\U\\2", formatted_helper_df_name, perl = TRUE)
  formatted_helper_df_name <- sub("(^|\\s)([a-z])", "\\1\\U\\2", formatted_helper_df_name, perl = TRUE)

  # Save the results to CSV files
  save_results(top_5_samples, helper_df_name, output_directory)

  # Save taxa data
  taxa_data <- taxa[taxa_asvs %in% unique(helper_df$ASV), ]
  write.csv(taxa_data, file.path(output_directory, paste0(helper_df_name, "_helper_taxa.csv")), row.names = TRUE)

  ## ---- Plot Top 5 Samples----
  ggplot(top_5_samples_long, aes(x = Sample, y = Abundance, fill = Helper_Label)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = "white"),
      plot.background = element_rect(fill = "white", color = "white"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = paste("Helper ASVs Abundance: Top 5 vs All Samples for", helper_df_name),
      x = "Samples", y = "Relative Abundance"
    ) +
    scale_x_discrete(drop = FALSE)

  # Save the plot
  ggsave(
    filename = paste0(output_directory, "/", helper_df_name, "_top_5_samples_abundance_plot.png"),
    width = 12, height = 6
  )

  ## ---- Plot All Samples----
  ggplot(all_samples_long, aes(x = Sample, y = Abundance, fill = Helper_Label)) +
    geom_bar(stat = "identity") +
    theme_classic() + # Ensures a clean white background
    theme(
      panel.grid = element_blank(), # Removes all grid lines
      panel.background = element_rect(fill = "white", color = "white"), # White background
      plot.background = element_rect(fill = "white", color = "white") # White outer background
    ) +
    labs(
      title = paste("All Samples for", formatted_helper_df_name, "Helper ASVs"),
      x = "Samples", y = "Relative Abundance"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(drop = FALSE)

  ggsave(
    filename = file.path(output_directory, paste0(helper_df_name, "_all_samples_abundance_plot.png")),
    width = 12, height = 6
  )

  cat("Processing complete.\n")
}
