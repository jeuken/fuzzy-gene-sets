import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

# Activate the automatic conversion between R and pandas DataFrame objects
pandas2ri.activate()

# Load R packages
ro.r('''
    library(tibble)
    library(limma)
    library(edgeR)
''')

# Set input and output paths
input_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/gemma/Alzheimers_GSE95587_formatted.csv"
output_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/gemma/Alzheimers_GSE95587_ranked.csv"

# Define R script to run limma analysis
r_script = f'''
    # Load necessary libraries
    library(tibble)
    library(limma)
    library(edgeR)

    # Load log2-transformed and quantile-normalized data
    data_matrix <- read.csv("{input_path}", row.names=1)

    # Extract column names
    column_names <- colnames(data_matrix)

    # Create a data frame with sample information
    sample_info <- data.frame(
      sample = column_names,
      condition = ifelse(grepl("^control", column_names), "control", "disease")
    )

    # Create the design matrix based on sample information
    design <- model.matrix(~ condition, data=sample_info)

    # Fit linear model to the log2-transformed data
    fit <- lmFit(data_matrix, design)

    # Apply empirical Bayes moderation
    fit <- eBayes(fit, trend = TRUE)

    # Extract results and adjust for multiple testing
    results <- topTable(fit, adjust="fdr", number=Inf)

    # Rename columns
    colnames(results)[colnames(results) == "P.Value"] <- "p_value"
    colnames(results)[colnames(results) == "adj.P.Val"] <- "adj_p_value"

    # Turn row names into Ensembl_ID column
    results <- rownames_to_column(results, "Ensembl_ID")

    # Sort results by descending t-statistic
    results_sorted <- results[order(results$t, decreasing=TRUE), ]

    # Save results to a CSV file
    write.csv(results_sorted, "{output_path}", row.names=FALSE)

    # Optional: View the top results
    head(results_sorted)
'''

# Execute the R script
ro.r(r_script)

print(f"Analysis complete. Results saved to {output_path}")
