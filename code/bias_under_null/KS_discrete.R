# Load necessary packages
library(readr)
library(dgof)

# Function to load p-values from a CSV file
load_p_values <- function(file_path) {
  # Load the file using read_csv for tab-separated values
  p_values_df <- read_csv(file_path, col_types = cols())
  return(p_values_df)
}

# Function to perform discrete KS test on all p-values
ks_test_all_p_values <- function(p_values) {
  if (length(p_values) > 0) {
    # Perform the discrete Kolmogorov-Smirnov test
    ks_test_result <- dgof::dks.test(p_values, punif)
    return(list(statistic = ks_test_result$statistic, p_value = ks_test_result$p.value))
  } else {
    return(list(statistic = NA, p_value = NA))  # Not enough data
  }
}

# Function to perform KS test per pathway
ks_test_per_pathway <- function(p_values_df) {
  results <- list()
  for (pathway in rownames(p_values_df)) {
    # Extract p-values for the pathway and remove NAs
    pathway_p_values <- na.omit(as.numeric(p_values_df[pathway, ]))
    
    # Only perform the test if there is at least one p-value
    if (length(pathway_p_values) > 0) {
      ks_result <- ks_test_all_p_values(pathway_p_values)
      results[[pathway]] <- list(
        KS_Statistic = ks_result$statistic, 
        p_value = ks_result$p_value
      )
    } else {
      results[[pathway]] <- list(KS_Statistic = NA, p_value = NA)
    }
  }
  return(results)
}

# Main function to determine uniformity
determine_uniformity <- function(p_values_file, output_file) {
  # Load p-values
  p_values_df <- load_p_values(p_values_file)
  
  # Convert all values to numeric, handling non-numeric values as NA
  p_values_df[] <- lapply(p_values_df, function(x) as.numeric(as.character(x)))
  
  # Perform KS test on all p-values
  all_p_values <- na.omit(unlist(p_values_df))
  ks_all_results <- ks_test_all_p_values(all_p_values)
  
  # Perform KS test per pathway
  ks_results_per_pathway <- ks_test_per_pathway(p_values_df)
  
  # Create summary data frame
  summary_df <- data.frame(
    Pathway_Name = names(ks_results_per_pathway),
    KS_Statistic = sapply(ks_results_per_pathway, function(x) x$KS_Statistic),
    p_value = sapply(ks_results_per_pathway, function(x) x$p_value)
  )
  
  # Add overall KS test results as an additional row
  summary_df <- rbind(summary_df, data.frame(Pathway_Name = "Overall", 
                                             KS_Statistic = ks_all_results$statistic, 
                                             p_value = ks_all_results$p_value))
  
  # Save the results to a CSV file
  write_csv(summary_df, output_file)
  print(paste("K-S test results saved to:", output_file))
}

# Define file paths
p_values_file <- '../../data/bias_under_null/null_distributions.csv'  # Input file with p-values
output_file <- '../../data/bias_under_null/ks_test_results.csv'  # Output file for K-S test results

# Run the main function
determine_uniformity(p_values_file, output_file)
