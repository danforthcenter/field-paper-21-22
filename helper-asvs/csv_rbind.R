library(dplyr)
library(readr)

# Set the parent directory
parent_dir <- "/Users/eflom/Library/CloudStorage/OneDrive-DDPSC/bart_lab/21_22_field_paper/Correlation/spearman_correlation/rr_samples"

# Get all CSV file paths from subdirectories
csv_files <- list.files(parent_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)

# Read all CSVs into a list of data frames
df_list <- lapply(csv_files, read_csv)

# Combine them into one data frame
combined_df <- bind_rows(df_list)

# Write the final CSV to the parent directory
write_csv(combined_df, file.path(parent_dir, "combined_output.csv"))

cat("CSV files successfully merged into combined_output.csv\n")
