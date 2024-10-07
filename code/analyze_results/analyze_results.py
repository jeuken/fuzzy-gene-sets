#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt

def load_ora_results(file_path):
    """Load ORA results from a CSV file."""
    return pd.read_csv(file_path)

def load_pathway_data(pathway_file):
    df = pd.read_csv(pathway_file, sep='\t', usecols=['Pathway_Name', 'Ensembl_ID', 'Overlap_Membership', 'Description'])
    
    # Check the loaded DataFrame
    print("Initial DataFrame columns:", df.columns)
    
    # Calculate the number of unique genes per pathway (Pathway_Size)
    pathway_size_df = df.groupby('Pathway_Name').agg(Pathway_Size=('Ensembl_ID', 'nunique')).reset_index()
    
    # Calculate the range of Overlap_Membership_linear for each pathway
    overlap_range_df = df.groupby('Pathway_Name')['Overlap_Membership'].agg(
        Overlap_Range=lambda x: x.max() - x.min()).reset_index()
    
    # Merge both results (pathway size and overlap range)
    pathway_data_df = pd.merge(pathway_size_df, overlap_range_df, on='Pathway_Name')

    # Add description
    description_df = df[['Pathway_Name', 'Description']].drop_duplicates()
    pathway_data_df = pd.merge(pathway_data_df, description_df, on='Pathway_Name')

    # Add Ensembl_ID as a column instead of an index
    # To do this, we need to ensure we include all relevant information
    ensembl_ids_df = df[['Pathway_Name', 'Ensembl_ID']].drop_duplicates()

    # Merge Ensembl IDs into the final DataFrame
    pathway_data_df = pd.merge(pathway_data_df, ensembl_ids_df, on='Pathway_Name', how='left')

    # Check the final DataFrame structure
    print("Pathway DataFrame columns after processing:", pathway_data_df.columns)  # Check columns again
    
    return pathway_data_df


def calculate_jaccard_index(pathway_data_df, reference_pathway):
    """Calculate Jaccard index for pathways compared to a reference pathway."""
    # Get genes from the reference pathway
    reference_genes = set(pathway_data_df[pathway_data_df['Pathway_Name'] == reference_pathway]['Ensembl_ID'].values)

    jaccard_results = []
    pathways_seen = set()  # To keep track of pathways already processed

    for _, row in pathway_data_df.iterrows():
        pathway_name = row['Pathway_Name']
        
        # Only process the pathway if it hasn't been seen before
        if pathway_name not in pathways_seen:
            pathways_seen.add(pathway_name)  # Mark this pathway as seen
            pathway_genes = set(pathway_data_df[pathway_data_df['Pathway_Name'] == pathway_name]['Ensembl_ID'].values)
            intersection = reference_genes.intersection(pathway_genes)
            union = reference_genes.union(pathway_genes)
            jaccard_index = len(intersection) / len(union) if union else 0
            
            jaccard_results.append({
                'Pathway_Name': pathway_name,
                'Jaccard_Index': jaccard_index,
                'Description': row['Description']
            })

    # Create a DataFrame and sort it
    jaccard_df = pd.DataFrame(jaccard_results)
    jaccard_df_sorted = jaccard_df.sort_values(by='Jaccard_Index', ascending=False).reset_index(drop=True)
    
    print(jaccard_df_sorted)  # Print or return sorted results as needed
    return jaccard_df_sorted



def plot_scatter_comparison(crisp_df, fuzzy_df, pathway_data_df):
    """Plot comparison of crisp and fuzzy ranks along with pathway characteristics."""
    merged_df = pd.merge(crisp_df[['Pathway_Name', 'Rank', 'p_value']], 
                           fuzzy_df[['Pathway_Name', 'Rank', 'p_value']], 
                           on='Pathway_Name', suffixes=('_Crisp', '_Fuzzy'))

    # Scatter plot: Crisp vs. Fuzzy Rank
    plt.figure(figsize=(8, 6))
    plt.scatter(merged_df['Rank_Crisp'], merged_df['Rank_Fuzzy'], alpha=0.6, color='b')
    plt.title('Crisp Rank vs Fuzzy Rank')
    plt.xlabel('Crisp Rank')
    plt.ylabel('Fuzzy Rank')
    plt.grid(True)
    plt.show()

    merged_df['Rank_Difference'] = merged_df['Rank_Fuzzy'] - merged_df['Rank_Crisp']
    merged_df = pd.merge(merged_df, pathway_data_df, on='Pathway_Name', how='left')

    plt.figure(figsize=(8, 6))
    plt.scatter(merged_df['Pathway_Size'], merged_df['Rank_Difference'], alpha=0.6, color='g')
    plt.title('Rank Difference vs Pathway Size')
    plt.xlabel('Pathway Size')
    plt.ylabel('Rank Difference (Crisp - Fuzzy)')
    plt.grid(True)
    plt.show()

    plt.figure(figsize=(8, 6))
    plt.scatter(merged_df['Overlap_Range'], merged_df['Rank_Difference'], alpha=0.6, color='m')
    plt.title('Rank Difference vs Overlap Membership Range')
    plt.xlabel('Overlap Membership Range')
    plt.ylabel('Rank Difference (Crisp - Fuzzy)')
    plt.grid(True)
    plt.show()

    plt.figure(figsize=(8, 6))
    plt.scatter(merged_df['p_value_Crisp'], merged_df['p_value_Fuzzy'], alpha=0.6, color='r')
    plt.title('Crisp vs Fuzzy P-Value Comparison')
    plt.xlabel('Crisp P-Value')
    plt.ylabel('Fuzzy P-Value')
    plt.xscale('log')
    plt.yscale('log')
    plt.plot([1e-10, 1], [1e-10, 1], 'r--')  # Reference line
    plt.grid(True)
    plt.show()

    merged_df['Absolute_Rank_Difference'] = merged_df['Rank_Difference'].abs()
    top_increases = merged_df.nlargest(10, 'Rank_Difference')[['Pathway_Name', 'Description', 'Rank_Difference']]
    top_decreases = merged_df.nsmallest(10, 'Rank_Difference')[['Pathway_Name', 'Description', 'Rank_Difference']]
    
    print("Top 10 Pathways with Increased Rank:")
    print(top_increases)
    
    print("\nTop 10 Pathways with Decreased Rank:")
    print(top_decreases)

def calculate_rank_change(crisp_results, fuzzy_results):
    """Calculate change in ranks between crisp and fuzzy results."""
    merged_ranks = pd.merge(crisp_results[['Pathway_Name', 'Rank']], 
                             fuzzy_results[['Pathway_Name', 'Rank']], 
                             on='Pathway_Name', suffixes=('_Crisp', '_Fuzzy'))

    # Calculate the change in rank
    merged_ranks['Rank_Change'] = merged_ranks['Rank_Fuzzy'] - merged_ranks['Rank_Crisp']
    
    return merged_ranks

def plot_jaccard_vs_rank_change(jaccard_df, rank_change_df):
    """Plot Jaccard index against rank change."""
    # Merge the Jaccard index DataFrame with the rank change DataFrame
    merged_df = pd.merge(jaccard_df, rank_change_df[['Pathway_Name', 'Rank_Change']], on='Pathway_Name', how='inner')

    # Create scatter plot
    plt.figure(figsize=(10, 6))
    plt.scatter(merged_df['Jaccard_Index'], merged_df['Rank_Change'], alpha=0.7, color='purple')
    plt.title('Jaccard Index vs Change in Ranks')
    plt.xlabel('Jaccard Index')
    plt.ylabel('Change in Ranks (Fuzzy - Crisp)')
    plt.axhline(0, color='gray', linestyle='--')  # Reference line for no change
    plt.grid(True)
    plt.show()

# Define file paths
crisp_file = '/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/ORA_output/infection/HIV/HIV_E-GEOD-111727/standard_ora_results_HIV_E-GEOD-111727_Crisp_Membership_Crisp_Membership.csv'
fuzzy_file = '/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/ORA_output/infection/HIV/HIV_E-GEOD-111727/ora_results_HIV_E-GEOD-111727_Crisp_Membership_Overlap_Membership.csv' 
pathway_file = "../../data/pathways/KEGG/KEGG_2022_membership.tsv"

# Load the results
crisp_results = load_ora_results(crisp_file)
fuzzy_results = load_ora_results(fuzzy_file)
pathway_data_df = load_pathway_data(pathway_file)

reference_pathway = "hsa05170"

# Calculate Jaccard index
jaccard_index_sorted = calculate_jaccard_index(pathway_data_df, reference_pathway)

# Calculate rank change
rank_change_df = calculate_rank_change(crisp_results, fuzzy_results)

# Plot Jaccard index against rank change
plot_jaccard_vs_rank_change(jaccard_index_sorted, rank_change_df)