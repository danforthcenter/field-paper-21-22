# Load necessary libraries
library(readr)
library(ggplot2)
library(tidyverse)
library(Hmisc)
library(compositions) # For CLR transformation

# Read ASV abundance table
asv_data <- read_csv("/Users/eflom/Library/CloudStorage/OneDrive-DDPSC/bart_lab/21_22_field_paper/Abundance/asvTable_Abnd_25_or_more.csv")

# Convert to numeric matrix (ensure it’s non-negative)
asv_data_numeric <- as.data.frame(lapply(asv_data, as.numeric))

# Step 1: Apply CLR transformation
clr_transformed <- clr(asv_data_numeric + 0.1) # Add 0.1 to avoid log(0) issues
if (is.list(clr_transformed)) {
  clr_transformed <- do.call(cbind, clr_transformed) # Convert list to matrix
}

# Step 1: Prepare the CLR transformed data
asv_data_matrix <- as.matrix(clr_transformed)

# Get the number of ASVs
n <- ncol(asv_data_matrix)

# Initialize an empty vector to store correlations and p-values
cor_values <- numeric(n)
p_values <- numeric(n)

# Define the ASV of interest
ASV_of_interest <- "ASV18093"
asv_of_interest_index <- which(colnames(asv_data_matrix) == ASV_of_interest)

# Loop through each ASV and compute Spearman correlation with the ASV of interest
for (i in 1:n) {
  # Skip if it's the same ASV
  if (i != asv_of_interest_index) {
    # Check if there are shared (non-NA) samples between ASV of interest and ASV i
    shared_samples <- asv_data_numeric[, asv_of_interest_index] > 0 & asv_data_numeric[, i] > 0

    if (any(shared_samples)) {
      correlation_result <- cor.test(asv_data_matrix[shared_samples, asv_of_interest_index], asv_data_matrix[shared_samples, i], method = "spearman", exact = FALSE)
      cor_values[i] <- correlation_result$estimate
      p_values[i] <- correlation_result$p.value
    }
  }
}


# Combine the correlation values and p-values into a data frame
cor_results <- data.frame(
  ASV = colnames(asv_data_matrix),
  Correlation = cor_values,
  P_Value = p_values
)

# View the results for the top correlations (p < 0.05)
significant_results <- cor_results[cor_results$P_Value < 0.05, ]
head(significant_results)

# Prepare the data for plotting significant results
plot_data <- significant_results

ggplot(significant_results, aes(x = 1:nrow(significant_results), y = Correlation)) +
  geom_point(aes(color = P_Value < 0.05), alpha = 0.6) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  scale_color_manual(values = c("gray", "blue")) + # Blue for significant, gray for non-significant
  theme_minimal() +
  labs(
    title = "Spearman Correlation of ASV18093 with Other ASVs",
    x = "ASVs",
    y = "Correlation Coefficient"
  ) +
  scale_x_continuous(breaks = seq(1, nrow(significant_results), by = 5000)) + # Set breaks at intervals of 5000
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
