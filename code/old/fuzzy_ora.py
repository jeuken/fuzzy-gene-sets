#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.stats import hypergeom


def load_memberships(query_file, query_membership_type):
    """Loads and returns the query set memberships as a dictionary, filtering out NaN values."""
    
    # Read the TSV file into a DataFrame, using 'Ensembl_ID' as the index column
    query_df = pd.read_csv(query_file, sep='\t', index_col='Ensembl_ID')
    
    
    # Filter out rows where the membership value is NaN
    query_df = query_df[query_df[query_membership_type].notna()]
    
    # Convert the specified membership type column to floats and create a dictionary from it
    query_set = query_df[query_membership_type].astype(float).to_dict()
    
    # Return the resulting dictionary
    return query_set

def load_pathways(pathway_file, pathway_membership_type):
    """Loads and groups pathways data, filtering out NaN values."""
    
    # Read the TSV file into a DataFrame with specified data types for columns
    pathways_df = pd.read_csv(pathway_file, sep='\t', dtype={pathway_membership_type: float, 'Ensembl_ID': str})
    
    # Filter out rows where the membership value is NaN
    pathways_df = pathways_df[pathways_df[pathway_membership_type].notna()]
    
    # Group by 'Pathway_Name' and aggregate information
    pathways_grouped = pathways_df.groupby('Pathway_Name').agg(
        Ensembl_ID=('Ensembl_ID', list),                # List of Ensembl IDs for each pathway
        Memberships=(pathway_membership_type, list),    # List of membership values for each pathway
        Description=('Description', 'first')            # Description of the pathway (assuming first entry is representative)
    ).reset_index()
    
    # Calculate the number of genes in each pathway
    pathways_grouped['Num_Genes'] = pathways_grouped['Ensembl_ID'].apply(len)
    
    # Return the resulting DataFrame
    return pathways_grouped

def compute_fuzzy_intersection(query_set, pathway_genes, pathway_memberships):
    """Calculates the fuzzy intersection size between query set and a pathway."""
    # Calculate the fuzzy intersection by summing the minimum of the query membership and pathway membership for each gene
    return sum(min(query_set.get(gene, 0), membership) for gene, membership in zip(pathway_genes, pathway_memberships))

def shuffle_memberships(query_set):
    """Shuffles the membership values of the query set."""
     # Shuffle the membership values randomly
    shuffled_values = np.random.permutation(list(query_set.values()))
    
    # Create a new dictionary with the original gene identifiers and the shuffled membership values
    return dict(zip(query_set.keys(), shuffled_values))

def compute_null_distributions(query_set, pathways_grouped, iterations=1000):
    """Computes the null distributions of fuzzy intersection sizes by shuffling query set memberships."""
    # Initialize a dictionary to store the null distributions for each pathway
    null_distributions = {pathway: [] for pathway in pathways_grouped['Pathway_Name']}
    
    # Iterate over the number of iterations to generate null distributions
    for iteration in range(iterations):
        # Shuffle the memberships of the query set
        shuffled_query_set = shuffle_memberships(query_set)
        
        # Compute fuzzy intersection sizes for each pathway with the shuffled query set
        for pathway, genes, memberships in zip(pathways_grouped['Pathway_Name'], pathways_grouped['Ensembl_ID'], pathways_grouped['Memberships']):
            intersection_size = compute_fuzzy_intersection(shuffled_query_set, genes, memberships)
            null_distributions[pathway].append(intersection_size)
        
        # Print progress every 1000 iterations
        if (iteration+1) % 1000 == 0:
            print(f"Iteration {iteration + 1}/{iterations}")
    
    # Return the computed null distributions
    return null_distributions

def calculate_p_values(observed_scores, null_distributions):
    """Calculates p-values based on the observed scores and null distributions."""
    # Initialize a dictionary to store the p-values for each pathway
    p_values = {}
    
    # Iterate over the observed scores for each pathway
    for pathway, observed_score in observed_scores.items():
        # Retrieve the null distribution for the current pathway
        null_dist = null_distributions[pathway]
        
        # Calculate the p-value as the proportion of null distribution scores that are greater than or equal to the observed score
        p_value = (sum(score >= observed_score for score in null_dist)) / len(null_dist)
        
        # Store the computed p-value in the dictionary
        p_values[pathway] = p_value
    
    # Return the dictionary with p-values
    return p_values

def plot_null_distribution(pathway, observed_score, null_distribution, p_value, plot_path, query_membership_type, pathway_membership_type, dataset_name):
    """Plots a histogram of the null distribution with the observed score and p-value."""
  # Create a new figure for the plot
    plt.figure(figsize=(8, 6))
    
    # Plot histogram of the null distribution
    plt.hist(null_distribution, bins = 30, alpha=0.7, color='gray', edgecolor='black')
    
    # Add a vertical line for the observed score
    plt.axvline(observed_score, color='red', linestyle='--', linewidth=2, label=f'Observed Score = {observed_score:.2f}')
    
    # Set the title of the plot with pathway and membership type information
    plt.title(f'Null Distribution for {pathway}\n(Query: {query_membership_type}, Pathway: {pathway_membership_type})')
    
    # Label the x and y axes
    plt.xlabel('Fuzzy Intersection Size')
    plt.ylabel('Frequency')
    
    # Add a legend to the plot
    plt.legend()
    
    # Annotate the plot with the p-value
    plt.annotate(f'P-value = {p_value:.4f}', xy=(0.7, 0.9), xycoords='axes fraction', fontsize=12, color='red')
    
    # Adjust layout to fit all elements
    plt.tight_layout()
    
    # Remove invalid characters and trailing underscores from the pathway name for the filename
    pathway_safe_name = pathway.replace(' ', '_').replace('/', '_').rstrip('_')
    
    # Save the plot as a PNG file in the specified output directory
    plt.savefig(f"{plot_path}/{dataset_name}_ora_null_{query_membership_type}_{pathway_membership_type}_{pathway_safe_name}.png")
    
    # Close the plot to free up memory
    plt.close()

def save_pathway_results(results_df, pathways_grouped, output_file, observed_scores, p_values):
    """Saves all pathways with their descriptions, number of genes, intersection sizes, and p-values, sorted by p-value."""
    # Create a DataFrame from the observed scores and p-values, including pathway descriptions and number of genes
    results_df = pd.DataFrame([
        {'Pathway': pathway,
         'Description': pathways_grouped[pathways_grouped['Pathway_Name'] == pathway]['Description'].values[0],
         'Num_Genes': pathways_grouped[pathways_grouped['Pathway_Name'] == pathway]['Num_Genes'].values[0],
         'Intersection_Size': observed_scores[pathway],
         'P_Value': p_values[pathway]}
        for pathway in observed_scores
    ])
    
    # Sort the DataFrame by p-value in ascending order
    results_df = results_df.sort_values(by='P_Value')
    
    # Save the results DataFrame to a CSV file
    results_df.to_csv(output_file, index=False)
    
    # Print a message indicating that the results have been saved
    print(f"All pathways results saved to {output_file}.")

def fuzzy_ora(query_file, pathway_file, specific_pathway=None, 
              query_membership_type='Membership_Score', pathway_membership_type='Default_Membership', 
              output_file='ora_enrichment_scores.csv', iterations=1000, plot_path=None):
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
        plot_path (str, optional): Directory to save the histogram plots of null distributions. Set to None to skip plotting.

    Returns:
        None
    """
    try:
        # Extract the dataset name from the query_file path by splitting the filename on underscores
        dataset_name = os.path.basename(query_file).split('_')[0] + "_" + os.path.basename(query_file).split('_')[1]
    
        # Load the memberships from the query file
        query_set = load_memberships(query_file, query_membership_type)
        
        # Load and group pathways data from the pathway file
        pathways_grouped = load_pathways(pathway_file, pathway_membership_type)
    
        # Filter the pathways to include only the specified pathway, if given
        if specific_pathway:
            pathways_grouped = pathways_grouped[pathways_grouped['Pathway_Name'] == specific_pathway]
    
        # Compute observed fuzzy intersection sizes for each pathway
        observed_scores = {
            pathway: compute_fuzzy_intersection(query_set, genes, memberships)
            for pathway, genes, memberships in zip(
                pathways_grouped['Pathway_Name'], 
                pathways_grouped['Ensembl_ID'], 
                pathways_grouped['Memberships']
            )
        }
    
        # Compute null distributions by shuffling memberships and calculating fuzzy intersections
        null_distributions = compute_null_distributions(query_set, pathways_grouped, iterations)
    
        # Calculate p-values based on observed scores and null distributions
        p_values = calculate_p_values(observed_scores, null_distributions)
    
        # Define the path for saving the results CSV file
        output_file = os.path.join(os.path.dirname(output_file), f"{dataset_name}_ora_scores_{query_membership_type}_{pathway_membership_type}.csv")
        
        # Save the results to a CSV file
        save_pathway_results(pd.DataFrame(), pathways_grouped, output_file, observed_scores, p_values)
    
        # Only proceed with plotting if plot_path is not None
        if plot_path:
            # Prepare the directory for saving plot images
            plot_plot_path = os.path.join(plot_path, f'{query_membership_type}_{pathway_membership_type}')
            os.makedirs(plot_plot_path, exist_ok=True)
    
            # Plot and save the null distributions for each pathway
            for pathway in observed_scores:
                plot_null_distribution(
                    pathway, 
                    observed_scores[pathway], 
                    null_distributions[pathway], 
                    p_values[pathway], 
                    plot_plot_path, 
                    query_membership_type, 
                    pathway_membership_type, 
                    dataset_name
                )
            # Print a message indicating where the results and plots have been saved
            print(f"Enrichment scores, p-values, and histograms saved to {output_file} and {plot_plot_path}.")
        else:
            # Print a message indicating where the results have been saved without plots
            print(f"Enrichment scores and p-values saved to {output_file}. No plots were generated.")
        
    except Exception as e:
        # Print the error message if an exception occurs
        print(f"Error: {e}")


# Application
query_file = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/gemma/Alzheimers_GSE95587_query_t.csv"
pathway_file = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/KEGG_2022_entrez_with_membership_overlap.tsv"
query_membership_type = "Crisp_Membership"  
pathway_membership_type = "Overlap_Membership"
output_file = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/ORA_output/test_ora_scores.csv"
plot_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/ORA_output"

# Run the main function
fuzzy_ora(query_file, pathway_file, query_membership_type=query_membership_type, 
          pathway_membership_type=pathway_membership_type, output_file=output_file, 
          iterations=1000, plot_path=plot_path)

