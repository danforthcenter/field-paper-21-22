#!/usr/bin/env Rscript
# Load necessary libraries
library(readr)
library(ggplot2)
library(tidyverse)
library(Hmisc)

# Get the ASV of interest from command-line arguments
args <- commandArgs(trailingOnly = TRUE)
ASV_of_interest <- args[1]
asv_data_file <- read_csv(args[2])
out_path <- args[3]

# Check the number of arguments
if (length(args) != 3) {
  stop(paste("Incorrect number of args"))
}

# Read ASV abundance table **once** at the beginning
asv_data_numeric <- as.data.frame(asv_data_file, show_col_types = TRUE)  # Check types
rownames(asv_data_numeric) <- asv_data_numeric[[1]]
asv_data_numeric <- asv_data_numeric[, -1]

# Debugging: Print first few rows of the data to inspect its structure
# cat("Inspecting the ASV data frame:\n")
# print(head(asv_data_numeric))  # Check the first few rows of the data

# Check if ASV_of_interest exists in columns
if (!(ASV_of_interest %in% colnames(asv_data_numeric))) {
  stop(paste("Error: ASV", ASV_of_interest, "not found in the dataset."))
}

# Extract samples for the ASV of interest
ASV_of_interest_samples <- asv_data_numeric[asv_data_numeric[[ASV_of_interest]] > 0, ]
cat("Number of rows in ASV_of_interest_samples:", nrow(ASV_of_interest_samples), "\n")


# Check if ASV_of_interest_samples is empty
if (nrow(ASV_of_interest_samples) == 0) {
  message(paste("No samples contain the ASV of interest:", ASV_of_interest), call. = FALSE)
  quit(status = 1, save = "no")
}

# Function to analyze ASV correlation and plot positive correlations
analyze_asv_correlation <- function(ASV_of_interest, base_save_path) {
  
  # Define save path
  save_path <- file.path(base_save_path, ASV_of_interest)
  
  # Check if directory exists
  if (dir.exists(save_path)) {
    message("Skipping ASV ", ASV_of_interest, ": Directory already exists - ", save_path)
    return(NULL)
  } else {
    # Create directory if it doesn't exist
    dir.create(save_path, recursive = TRUE)
    message("Created new directory: ", save_path)
  }
  
  # Identify ASVs that appear in at least one of these samples (union of nonzero ASVs)
  union_ASVs <- colnames(ASV_of_interest_samples)[colSums(ASV_of_interest_samples) > 0]
  
  # Subset ASV_of_interest_shared_ASVs to include only ASVs in the union
  ASV_union_subset <- ASV_of_interest_samples[, union_ASVs]
  
  asv_of_interest_index <- which(colnames(ASV_union_subset) == ASV_of_interest)
  
  # Get the number of ASVs
  n <- ncol(ASV_union_subset)
  
  # Initialize result vectors
  cor_values <- rep(NA, n)
  names(cor_values) <- colnames(ASV_union_subset)
  
  p_values <- rep(NA, n)
  names(p_values) <- colnames(ASV_union_subset)
  
  # Loop through all other ASVs
  for (i in 1:n) {
    if (i != asv_of_interest_index) {
      # Check if ASV has nonzero values
      if (any(ASV_union_subset[, i]) > 0) {
        # Compute Spearman correlation
        correlation_result <- cor.test(ASV_union_subset[, asv_of_interest_index], 
                                       ASV_union_subset[, i], 
                                       method = "spearman", exact = FALSE)
        cor_values[i] <- correlation_result$estimate
        p_values[i] <- correlation_result$p.value
      }
    }
  }
  
  # Combine results into a data frame
  cor_results <- data.frame(ASV = colnames(ASV_union_subset),
                            Correlation = cor_values,
                            P_Value = p_values)
  
  # Filter for significant correlations (p < 0.05)
  significant_results <- cor_results[cor_results$P_Value < 0.05, ]
  
  # Save correlation results as a CSV file in the ASV directory
  tryCatch({
    write_csv(significant_results, file.path(save_path, paste0("correlation_results_", ASV_of_interest, ".csv")))
  }, error = function(e) {
    message("Error saving correlation results: ", e)
  })
  
  # Filter for positive correlations
  positive_results <- significant_results[significant_results$Correlation > 0, ]
  
  # Loop through each positive correlation and plot
  for (i in 1:nrow(positive_results)) {
    # Get the name of the ASV with positive correlation
    asv_name <- positive_results$ASV[i]
    
    # Extract the Spearman correlation coefficient
    spearman_corr <- positive_results$Correlation[i]
    
    # Check if both ASV_of_interest and correlated ASV exist as columns
    if (asv_name %in% colnames(ASV_union_subset) && ASV_of_interest %in% colnames(ASV_union_subset)) {
      # Extract the columns for ASV of interest and the correlated ASV
      ASV_of_interest_column <- ASV_union_subset[[ASV_of_interest]]
      correlated_asv_column <- ASV_union_subset[[asv_name]]
      
      # Check for NA values in the columns
      if (any(is.na(ASV_of_interest_column)) || any(is.na(correlated_asv_column))) {
        warning(paste("NA values found for", ASV_of_interest, "or", asv_name))
        next  # Skip this iteration if there are NA values
      }
      
      # Create the data frame for plotting
      plot_data <- data.frame(
        asvi = ASV_of_interest_column,
        Correlated_ASV = correlated_asv_column
      )
      
      # Fit a linear model for the data
      lm_model <- lm(Correlated_ASV ~ asvi, data = plot_data)
      
      # Extract the R² value from the model
      r_squared <- summary(lm_model)$r.squared
      
      # Create scatter plot with line of best fit
      p <- ggplot(plot_data, aes(x = asvi, y = Correlated_ASV)) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add line of best fit
        labs(title = paste("Correlation between", ASV_of_interest, "and", asv_name),
             x = ASV_of_interest,
             y = asv_name) +
        annotate("text", x = max(plot_data$asvi) * 0.8, y = max(plot_data$Correlated_ASV) * 0.9,
                 label = paste("Spearman's rho = ", round(spearman_corr, 2)), color = "red", size = 5) +  # Add Spearman's rho on the plot
        theme_minimal()
      
      # Save plot
      file_name <- sprintf("%s/scatter_plot_%s_vs_%s.png", save_path, ASV_of_interest, asv_name)
      ggsave(file_name, plot = p, width = 8, height = 6, dpi = 300)
    }
  }
  
  # Return results
  #return(list(ASV_union_subset = ASV_union_subset, Significant_Correlations = significant_results, Save_Path = save_path))
}

# Call the function to analyze ASV correlation
analyze_asv_correlation(ASV_of_interest, out_path)
