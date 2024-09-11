#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from rpy2.robjects import pandas2ri, r
from rpy2.robjects.packages import importr
import rpy2.robjects as ro

# Activate the conversion from pandas to R dataframe
pandas2ri.activate()

# Import R's base and limma libraries
base = importr('base')
limma = importr('limma')

def shuffle_phenotype_labels(df):
    """Randomly shuffles phenotype labels among 'control' and 'disease' columns."""
    shuffled_df = df.copy()
    
    # Ensure 'Ensembl_ID' is a column
    if shuffled_df.index.name == 'Ensembl_ID':
        shuffled_df.reset_index(inplace=True)
    
    phenotype_cols = [col for col in shuffled_df.columns if col.startswith(('control', 'disease'))]
    shuffled_col_names = np.random.permutation(phenotype_cols)
    shuffled_df = shuffled_df.rename(columns=dict(zip(phenotype_cols, shuffled_col_names)))
    
    return shuffled_df

def calculate_statistics_from_expression_dataframe(df):
    """Calculates statistics from expression dataframe using limma in R."""
    # Check for duplicates and missing values
    duplicates = df['Ensembl_ID'].duplicated().sum()
    print(f"Number of duplicate Ensembl_IDs: {duplicates}")
    print(f"Number of missing Ensembl_IDs: {df['Ensembl_ID'].isnull().sum()}")
    
    # Set Ensembl_ID as the row index
    df.set_index('Ensembl_ID', inplace=True)
    
    # Convert the pandas DataFrame to an R DataFrame
    r_df = pandas2ri.py2rpy(df)

    # Assign the R DataFrame to the R environment
    ro.globalenv['data_matrix'] = r_df

    # Define the R script to run limma analysis
    r_script = '''
        # Load necessary libraries
        library(tibble)
        library(limma)
        library(edgeR)
    
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
        
        return(results_sorted)
    '''
    
    # Execute the R script
    ranked_results = ro.r(r_script)
    
    # Convert the R dataframe back to a pandas dataframe
    results_df = pandas2ri.rpy2py(ranked_results)

    # Ensure Ensembl_ID is a column
    if results_df.index.name == 'Ensembl_ID':
        results_df.reset_index(inplace=True)

    return results_df

def fuzzy_gsea_score(ranked_list, pathways_dict):
    """Calculates enrichment scores for pathways using a ranked list and a dataset of pathways."""
    enrichment_scores = {}
    plotting_values_all = {}
    N = len(ranked_list)
    ranked_gene_list = list(ranked_list.keys())

    for pathway_name, pathway_info in pathways_dict.items():
        gene_ids = set(pathway_info['genes'].keys())
        correlation = {gene_id: ranked_list.get(gene_id, 0) for gene_id in gene_ids}
        N_r = sum(abs(correlation.get(gene_id, 0)) * pathway_info['genes'].get(gene_id, 1) for gene_id in gene_ids)
        N_misses = N - len(gene_ids)
        P_hit = 0
        P_miss = 1
        counter = 0
        enrichment_score = 0.0
        plotting_values = []

        for idx, gene_id in enumerate(ranked_gene_list):
            if gene_id in gene_ids:
                membership_value = pathway_info['genes'].get(gene_id, 1)
                P_hit += abs(correlation.get(gene_id, 0) * membership_value) / N_r
                counter += 1
            P_miss = ((idx - counter) + 1) / N_misses

            # Update enrichment score if the current score is higher
            if abs(P_hit - P_miss) > abs(enrichment_score):
                enrichment_score = P_hit - P_miss

            # Track the enrichment score for plotting
            plotting_values.append(P_hit - P_miss)

        enrichment_scores[pathway_name] = enrichment_score
        plotting_values_all[pathway_name] = plotting_values

    return enrichment_scores, plotting_values_all

def permute_and_calculate_null_distribution(expression_df, pathways_dict, n_permutations=1000):
    """Generates a null distribution of enrichment scores by permuting column labels."""
    null_distributions = {pathway: [] for pathway in pathways_dict.keys()}
    
    for i in range(n_permutations):
        shuffled_df = shuffle_phenotype_labels(expression_df)
        ranked_df = calculate_statistics_from_expression_dataframe(shuffled_df)
        ranked_list = ranked_df.set_index('Ensembl_ID')['t'].to_dict()
        enrichment_scores, _ = fuzzy_gsea_score(ranked_list, pathways_dict)

        for pathway, score in enrichment_scores.items():
            null_distributions[pathway].append(score)
        
        if (i + 1) % 100 == 0:
            print(f'Completed {i + 1} permutations.')
    
    return null_distributions

def calculate_p_value(observed_score, null_distribution):
    """Calculates the two-sided p-value for an observed score given a null distribution."""
    count = sum(1 for score in null_distribution if abs(score) >= abs(observed_score))
    p_value = count / len(null_distribution)
    return p_value

def save_enrichment_plots(plotting_values_all, observed_scores, pathways_dict, plot_path, membership):
    """Saves enrichment plots for each pathway with titles including the pathway description."""
    membership_path = os.path.join(plot_path, membership, 'enrichment_plots')
    os.makedirs(membership_path, exist_ok=True)

    for pathway, plotting_values in plotting_values_all.items():
        plt.figure(figsize=(10, 6))
        plt.plot(range(len(plotting_values)), plotting_values, label='Enrichment Score')
        final_score = observed_scores[pathway]
        plt.axhline(y=final_score, color='r', linestyle='--', label=f'Final Enrichment Score: {final_score:.2f}')
        pathway_description = pathways_dict[pathway]['description']
        plt.title(f'GSEA Plot for {pathway}\nDescription: {pathway_description}\nMembership: {membership}')
        plt.xlabel('Rank')
        plt.ylabel('Enrichment Score')
        plt.legend()
        plt.grid(True)
        plt.savefig(os.path.join(membership_path, f'{pathway}_{membership}_gsea_plot.png'))
        plt.close()

def save_null_distribution_plots(null_distributions, observed_scores, pathways_dict, plot_path, membership):
    """Saves histograms of the null distributions with the observed enrichment score, p-value, and pathway description."""
    null_dist_path = os.path.join(plot_path, membership, 'null_distributions')
    os.makedirs(null_dist_path, exist_ok=True)

    for pathway, null_distribution in null_distributions.items():
        observed_score = observed_scores[pathway]
        p_value = calculate_p_value(observed_score, null_distribution)
        pathway_description = pathways_dict[pathway]['description']
        plt.figure(figsize=(10, 6))
        plt.hist(null_distribution, bins=30, alpha=0.75, color='blue', label='Null Distribution')
        plt.axvline(observed_score, color='red', linestyle='--', linewidth=2, label=f'Observed Score: {observed_score:.2f}')
        plt.title(f'Null Distribution for {pathway}\nDescription: {pathway_description}\nMembership: {membership}')
        plt.xlabel('Enrichment Score')
        plt.ylabel('Frequency')
        plt.text(0.95, 0.9, f'p-value: {p_value:.4f}', transform=plt.gca().transAxes, 
                 fontsize=12, verticalalignment='top', horizontalalignment='right')
        plt.legend()
        plt.grid(True)
        plt.savefig(os.path.join(null_dist_path, f'{pathway}_{membership}_null_distribution.png'))
        plt.close()

def main(expression_path, pathway_file_path, output_path, plot_path=None, membership='Default_Membership', n_permutations=1000):
    expression_df = pd.read_csv(expression_path)                     
    ranked_df = calculate_statistics_from_expression_dataframe(expression_df)
    ranked_list = ranked_df.set_index('Ensembl_ID')['t'].to_dict()
    
    pathways_df = pd.read_csv(pathway_file_path, sep='\t')
    pathways_df = pathways_df.dropna(subset=[membership])
    pathways_dict = {}

    for _, row in pathways_df.iterrows():
        pathway_name = row['Pathway_Name']
        gene_id = row['Ensembl_ID']
        membership_value = row[membership]
        pathway_description = row.get('Description', '')
        if pathway_name not in pathways_dict:
            pathways_dict[pathway_name] = {'genes': {}, 'description': pathway_description}
        pathways_dict[pathway_name]['genes'][gene_id] = membership_value

    observed_scores, plotting_values_all = fuzzy_gsea_score(ranked_list, pathways_dict)
    null_distributions = permute_and_calculate_null_distribution(expression_df, pathways_dict, n_permutations)
    
    results = []

    for pathway, observed_score in observed_scores.items():
        null_distribution = null_distributions[pathway]
        p_value = calculate_p_value(observed_score, null_distribution)
        pathway_description = pathways_dict[pathway]['description']
        results.append({'Pathway_Name': pathway, 'Observed_Score': observed_score, 'p-value': p_value, 'Description': pathway_description})

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values(by='p-value')  # Rank results by p-value
    results_df.to_csv(os.path.join(output_path, f'fuzzy_gsea_results_{membership}.tsv'), sep='\t', index=False)

    if plot_path:
        save_enrichment_plots(plotting_values_all, observed_scores, pathways_dict, plot_path, membership)
        save_null_distribution_plots(null_distributions, observed_scores, pathways_dict, plot_path, membership)

if __name__ == "__main__":
    expression_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/gemma/Alzheimers_GSE95587_formatted.csv"
    pathway_file_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/KEGG_2022_entrez_with_membership.tsv"
    output_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/GSEA_output"
    plot_path = '/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/GSEA_output/plots'
    main(expression_path, pathway_file_path, output_path, plot_path, membership='Overlap_Membership', n_permutations=10)
