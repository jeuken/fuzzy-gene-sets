#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def load_memberships(query_file, query_membership_type):
    """Loads and returns the query set memberships as a dictionary."""
    query_set = pd.read_csv(query_file, sep='\t', index_col='Ensembl_ID')[query_membership_type].astype(float).to_dict()
    return query_set

def load_pathways(pathway_file, pathway_membership_type):
    """Loads and groups pathways data."""
    pathways_df = pd.read_csv(pathway_file, sep='\t', dtype={pathway_membership_type: float, 'Ensembl_ID': str})
    pathways_grouped = pathways_df.groupby('Pathway_Name').agg(
        Ensembl_ID=('Ensembl_ID', list), 
        Memberships=(pathway_membership_type, list),
        Description=('Description', 'first')
    ).reset_index()
    return pathways_grouped

def compute_fuzzy_intersection(query_set, pathway_genes, pathway_memberships):
    """Calculates the (fuzzy) intersection size between query set and a pathway."""
    return sum(min(query_set.get(gene, 0), membership) for gene, membership in zip(pathway_genes, pathway_memberships))

def shuffle_memberships(query_set):
    """Shuffles the membership values of the query set."""
    shuffled_values = np.random.permutation(list(query_set.values()))
    return dict(zip(query_set.keys(), shuffled_values))

def compute_null_distributions(query_set, pathways_grouped, iterations=1000):
    """Computes the null distributions of fuzzy intersection sizes by shuffling query set memberships."""
    null_distributions = {pathway: [] for pathway in pathways_grouped['Pathway_Name']}
    
    for iteration in range(iterations):
        shuffled_query_set = shuffle_memberships(query_set)
        for pathway, genes, memberships in zip(pathways_grouped['Pathway_Name'], pathways_grouped['Ensembl_ID'], pathways_grouped['Memberships']):
            intersection_size = compute_fuzzy_intersection(shuffled_query_set, genes, memberships)
            null_distributions[pathway].append(intersection_size)
        
        # Print the number of pathways processed at each iteration
        if iteration % 1000 == 0:
            print(f"Iteration {iteration + 1}/{iterations}")
    
    return null_distributions

def calculate_p_values(observed_scores, null_distributions):
    """Calculates p-values based on the observed scores and null distributions."""
    p_values = {}
    for pathway, observed_score in observed_scores.items():
        null_dist = null_distributions[pathway]
        p_value = (sum(score >= observed_score for score in null_dist)) / (len(null_dist))
        p_values[pathway] = p_value
    return p_values

def plot_null_distribution(pathway, observed_score, null_distribution, p_value, output_dir, query_membership_type, pathway_membership_type, dataset_name):
    """Plots a histogram of the null distribution with the observed score and p-value."""
    plt.figure(figsize=(8, 6))
    plt.hist(null_distribution, bins=30, alpha=0.7, color='gray', edgecolor='black')
    plt.axvline(observed_score, color='red', linestyle='--', linewidth=2, label=f'Observed Score = {observed_score:.2f}')
    plt.title(f'Null Distribution for {pathway}\n(Query: {query_membership_type}, Pathway: {pathway_membership_type})')
    plt.xlabel('Fuzzy Intersection Size')
    plt.ylabel('Frequency')
    plt.legend()
    plt.annotate(f'P-value = {p_value:.4f}', xy=(0.7, 0.9), xycoords='axes fraction', fontsize=12, color='red')
    plt.tight_layout()
    
    # Remove trailing underscore from pathway name
    pathway_safe_name = pathway.replace(' ', '_').replace('/', '_').rstrip('_')
    plt.savefig(f"{output_dir}/{dataset_name}_ora_null_{query_membership_type}_{pathway_membership_type}_{pathway_safe_name}.png")
    plt.close()

def save_pathway_results(results_df, pathways_grouped, output_file, observed_scores, p_values):
    """Saves all pathways with their descriptions, intersection sizes, and p-values, sorted by p-value."""
    results_df = pd.DataFrame([
        {'Pathway': pathway, 'Description': pathways_grouped[pathways_grouped['Pathway_Name'] == pathway]['Description'].values[0],
         'Intersection_Size': observed_scores[pathway], 'P_Value': p_values[pathway]}
        for pathway in observed_scores
    ])
    results_df = results_df.sort_values(by='P_Value')
    results_df.to_csv(output_file, index=False)  # Save as CSV file
    print(f"All pathways results saved to {output_file}.")

def fuzzy_ora(query_file, pathway_file, specific_pathway=None, 
              query_membership_type='Membership_Score', pathway_membership_type='Default_Membership', 
              output_file='ora_enrichment_scores.csv', iterations=1000, output_dir='null_distributions'):
    """
    Main function to perform fuzzy ORA, compute null distributions, calculate p-values, 
    and save pathway results.

    Args:
        query_file (str): Path to the query set membership file.
        pathway_file (str): Path to the pathway information file.
        specific_pathway (str, optional): Specific pathway to focus on, or None for all pathways.
        query_membership_type (str, optional): Column name for query set membership.
        pathway_membership_type (str, optional): Column name for pathway membership.
        output_file (str, optional): Path to save the enrichment scores and p-values.
        iterations (int, optional): Number of shuffling iterations for null distribution computation.
        output_dir (str, optional): Directory to save the histogram plots of null distributions.

    Returns:
        None
    """
    try:
        # Extract the dataset name from the query_file path
        dataset_name = os.path.basename(query_file).split('_')[0] + "_" + os.path.basename(query_file).split('_')[1]

        # Ensure the output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # Load data
        query_set = load_memberships(query_file, query_membership_type)
        pathways_grouped = load_pathways(pathway_file, pathway_membership_type)

        # Filter to a specific pathway if specified
        if specific_pathway:
            pathways_grouped = pathways_grouped[pathways_grouped['Pathway_Name'] == specific_pathway]

        # Compute observed fuzzy intersection sizes
        observed_scores = {
            pathway: compute_fuzzy_intersection(query_set, genes, memberships)
            for pathway, genes, memberships in zip(
                pathways_grouped['Pathway_Name'], pathways_grouped['Ensembl_ID'], pathways_grouped['Memberships']
            )
        }

        # Compute null distributions
        null_distributions = compute_null_distributions(query_set, pathways_grouped, iterations)

        # Calculate p-values
        p_values = calculate_p_values(observed_scores, null_distributions)

        # Save results
        output_file = os.path.join(os.path.dirname(output_file), f"{dataset_name}_ora_scores_{query_membership_type}_{pathway_membership_type}.csv")
        save_pathway_results(pd.DataFrame(), pathways_grouped, output_file, observed_scores, p_values)

        # Prepare output directory for plots
        plot_output_dir = os.path.join(output_dir, f'{query_membership_type}_{pathway_membership_type}')
        os.makedirs(plot_output_dir, exist_ok=True)

        # Plot null distributions
        for pathway in observed_scores:
            plot_null_distribution(pathway, observed_scores[pathway], null_distributions[pathway], 
                                   p_values[pathway], plot_output_dir, query_membership_type, pathway_membership_type, dataset_name)

        print(f"Enrichment scores, p-values, and histograms saved to {output_file} and {plot_output_dir}.")

    except Exception as e:
        print(f"Error: {e}")

# Application
query_file = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/gemma/Alzheimers_GSE95587_query.csv"
pathway_file = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/KEGG_2022_entrez_with_membership.tsv"
query_membership_type = "Crisp_Membership"  
pathway_membership_type = "Expansion_Membership"

output_file = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/ORA_output/Alzheimers_GSE95587/ora_scores.csv"
output_dir = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/ORA_output/Alzheimers_GSE95587/null_distributions"

# Run the main function
fuzzy_ora(query_file, pathway_file, query_membership_type=query_membership_type, 
          pathway_membership_type=pathway_membership_type, output_file=output_file, 
          iterations=10000, output_dir=output_dir)

