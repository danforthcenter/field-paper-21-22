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
