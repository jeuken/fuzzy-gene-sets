#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def process_and_expand_pathways(gmt_path, Ensembl_ID_dict_path, string_scores_path, output_file):
    """
    Processes pathways from a GMT file, maps gene names to gene IDs using a provided gene ID dictionary,
    expands pathways with STRING scores, and filters and normalizes the expanded pathways.

    Args:
        gmt_path (str): Path to the GMT file containing pathway data.
        Ensembl_ID_dict_path (str): Path to the gene ID dictionary file.
        string_scores_path (str): Path to the STRING scores file.
        output_file (str): Path to save the final filtered and normalized pathways data.

    Returns:
        pd.DataFrame: Final DataFrame containing processed, expanded, and normalized pathways information.
    """
    try:
        # Load gene ID dictionary
        gene_dict = pd.read_csv(
            Ensembl_ID_dict_path,
            sep='\t',
            header=1,
            names=['Ensembl_ID', 'Gene_Name'],
            dtype={'Gene_Name': str}
        ).assign(Gene_Name=lambda df: df['Gene_Name'].str.strip()) \
         .dropna(subset=['Gene_Name']) \
         .drop_duplicates(subset='Gene_Name') \
         .set_index('Gene_Name')

        # Initialize lists to collect pathway data
        pathway_rows = []
        not_found_genes = set()

        # Read and process the GMT file
        with open(gmt_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                pathway_name, description, *gene_names = parts
                gene_names = [name.strip() for name in gene_names]
                for name in gene_names:
                    if name in gene_dict.index:
                        Ensembl_ID = gene_dict.loc[name, 'Ensembl_ID']
                        pathway_rows.append([
                            pathway_name,
                            description,
                            Ensembl_ID,
                            1  # Default membership
                        ])
                    else:
                        not_found_genes.add(name)

        print(f"Total number of pathway genes not found in dictionary: {len(not_found_genes)}")

        # Create DataFrame from collected rows
        pathways_df = pd.DataFrame(
            pathway_rows, 
            columns=['Pathway_Name', 'Description', 'Ensembl_ID', 'Crisp_Membership']
        )

        # Calculate gene frequencies for overlap membership
        gene_counts = pathways_df['Ensembl_ID'].value_counts()
        min_freq = gene_counts.min()
        max_freq = gene_counts.max()
        num_pathways = pathways_df['Pathway_Name'].nunique()

        print(f"Number of unique pathways: {num_pathways}")
        print(f"Maximum gene frequency: {max_freq}")
        print(f"Minimum gene frequency: {min_freq}")

        # Calculate overlap membership
        pathways_df['Overlap_Membership'] = (max_freq - pathways_df['Ensembl_ID'].map(gene_counts)) / (max_freq - min_freq)

        # Load and process STRING scores
        string_scores = pd.read_csv(string_scores_path, sep="\t")
        pathways_df['Expansion_Membership'] = 1

        # Add new genes from STRING scores
        new_rows = []
        completed_pathways_count = 0
        total_pathways = pathways_df['Pathway_Name'].nunique()

        for pathway_name, pathway_df in pathways_df.groupby('Pathway_Name'):
            pathway_genes = set(pathway_df['Ensembl_ID'])
            pathway_description = pathway_df['Description'].iloc[0]
            pathway_string = string_scores[string_scores['Gene2'].isin(pathway_genes)]
            gene1_scores = pathway_string.groupby('Gene1')['string_score'].mean()
            for gene, avg_score in gene1_scores.items():
                if gene not in pathway_genes:
                    new_row = {
                        'Pathway_Name': pathway_name,
                        'Ensembl_ID': gene,
                        'Expansion_Membership': avg_score / 1000,
                        'Description': pathway_description,
                        'Crisp_Membership': np.nan,  # Fill NaN for the new row
                        'Overlap_Membership': np.nan  # Fill NaN for the new row
                    }
                    new_rows.append(new_row)
            completed_pathways_count += 1
            print(f"Number of pathways completed: {completed_pathways_count}/{total_pathways}")

        # Create and save the final DataFrame
        new_rows_df = pd.DataFrame(new_rows)
        pathways_expanded = pd.concat([pathways_df, new_rows_df], ignore_index=True)
        pathways_expanded = pathways_expanded[pathways_expanded['Expansion_Membership'] >= 0.4]

        # Normalize Expansion_Membership
        min_value, max_value = 0.4, 1
        pathways_expanded['Expansion_Membership'] = pathways_expanded['Expansion_Membership'].apply(
            lambda x: (x - min_value) / (max_value - min_value) if pd.notna(x) else x
        )

        # Fill NaNs with zeros for memberships from original pathways data
        final_df = pathways_expanded
        # Save the final DataFrame
        final_df.to_csv(output_file, sep="\t", index=False)
        print(f"Filtered and normalized DataFrame saved to {output_file}")

        # Calculate and print max and min values below 1 in 'Expansion_Membership'
        filtered_df = final_df[final_df['Expansion_Membership'] < 1]
        print(f"Maximum value (below 1) in 'Expansion_Membership' : {filtered_df['Expansion_Membership'].max()}")
        print(f"Minimum value in 'Expansion_Membership': {filtered_df['Expansion_Membership'].min()}")

        # Create histograms for each pathway and save them
        for pathway_name, pathway_df in final_df.groupby('Pathway_Name'):
            # Ensure directories exist
            exp_dir = os.path.join(output_file.rsplit('/', 1)[0], 'Expansion_Membership')
            if not os.path.exists(exp_dir):
                os.makedirs(exp_dir)
            
            overlap_dir = os.path.join(output_file.rsplit('/', 1)[0], 'Overlap_Membership')
            if not os.path.exists(overlap_dir):
                os.makedirs(overlap_dir)
            
            # Plot and save Expansion_Membership histogram
            plt.figure(figsize=(10, 6))
            plt.hist(pathway_df['Expansion_Membership'].dropna(), bins=30, edgecolor='black', density=False)
            plt.title(f'Histogram of Expansion_Membership for {pathway_name}')
            plt.xlabel('Expansion_Membership')
            plt.ylabel('Frequency')
            plt.savefig(os.path.join(exp_dir, f'{pathway_name}_histogram.png'))
            plt.close()

            # Plot and save Overlap_Membership histogram
            plt.figure(figsize=(10, 6))
            # Exclude NaN values for histogram
            plt.hist(pathway_df['Overlap_Membership'].dropna(), bins=30, edgecolor='black', density=False)
            plt.title(f'Histogram of Overlap_Membership for {pathway_name}')
            plt.xlabel('Overlap_Membership')
            plt.ylabel('Frequency')
            plt.savefig(os.path.join(overlap_dir, f'{pathway_name}_histogram.png'))
            plt.close()

            print(f"Histograms for {pathway_name} saved to {exp_dir} and {overlap_dir}")

        return final_df

    except pd.errors.ParserError as e:
        print(f"Parsing Error: {e}")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

# Example usage
gmt_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/KEGG_2022_entrez.gmt"
Ensembl_ID_dict_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/IDs/entrez2ensembl.txt"
string_scores_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/STRING/STRING_GGI.txt"
output_file = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/KEGG_2022_entrez_with_membership.tsv"

process_and_expand_pathways(gmt_path, Ensembl_ID_dict_path, string_scores_path, output_file)
