#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

def process_and_save_pathways(gmt_path, gene_id_dict_path, output_path):
    """
    Processes pathways from a GMT file and maps gene names to gene IDs using a provided gene ID dictionary.
    Calculates default, random, and overlap memberships, and saves the results to a TSV file.

    Args:
        gmt_path (str): Path to the GMT file containing pathway data.
        gene_id_dict_path (str): Path to the gene ID dictionary file.
        output_path (str): Path to save the processed pathways data.

    Returns:
        pd.DataFrame: DataFrame containing processed pathway information.
    """
    try:
        # Load and preprocess the gene ID dictionary
        gene_dict = (
            pd.read_csv(
                gene_id_dict_path,
                sep='\t',
                header=1,
                names=['Ensembl_ID', 'Gene_Name'],
                dtype={'Gene_Name': str}
            )
            .assign(Gene_Name=lambda df: df['Gene_Name'].str.strip())  # Strip whitespace from gene names
            .dropna(subset=['Gene_Name'])  # Remove rows with missing gene names
            .drop_duplicates(subset='Gene_Name')  # Remove duplicate gene names
            .set_index('Gene_Name')  # Set gene names as index for mapping
        )

        # Initialize a list to collect pathway data
        rows = []
        not_found_genes = set()  # Set to track genes not found in the dictionary

        # Read and process the GMT file
        with open(gmt_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue  # Skip lines with insufficient data

                pathway_name, description, *gene_names = parts
                gene_names = [name.strip() for name in gene_names]  # Clean gene names

                # Map gene names to Ensembl IDs and assign memberships
                for name in gene_names:
                    if name in gene_dict.index:
                        gene_id = gene_dict.loc[name, 'Ensembl_ID']

                        # Assign membership values
                        default_membership = 1
                        random_membership = np.random.uniform(0, 1)  # Random membership between 0 and 1

                        # Append the row to the list
                        rows.append([
                            pathway_name,
                            description,
                            gene_id,
                            default_membership,
                            random_membership
                        ])
                    else:
                        not_found_genes.add(name)  # Track genes not found in the dictionary

        # Report the total number of genes not found
        print(f"Total number of pathway genes not found in dictionary: {len(not_found_genes)}")

        # Create DataFrame from collected rows
        pathways_df = pd.DataFrame(
            rows, 
            columns=['Pathway_Name', 'Description', 'Ensembl_ID', 'Default_Membership', 'Random_Membership']
        )

        # Calculate the frequency of each gene across pathways
        gene_counts = pathways_df['Ensembl_ID'].value_counts()

        # Get min and max frequencies for overlap membership calculation
        min_freq = gene_counts.min()
        max_freq = gene_counts.max()
        num_pathways = pathways_df['Pathway_Name'].nunique()

        # Print summary statistics
        print(f"Number of unique pathways: {num_pathways}")
        print(f"Maximum gene frequency: {max_freq}")
        print(f"Minimum gene frequency: {min_freq}")

        # Calculate overlap membership
        pathways_df['Gene_Frequency'] = pathways_df['Ensembl_ID'].map(gene_counts)
        pathways_df['Overlap_Membership'] = (max_freq - pathways_df['Gene_Frequency']) / (max_freq - min_freq)

        # Drop the temporary 'Gene_Frequency' column
        pathways_df.drop(columns=['Gene_Frequency'], inplace=True)

        # Save the processed DataFrame to a TSV file
        pathways_df.to_csv(output_path, sep='\t', index=False)
        print(f"Processed pathways saved successfully to {output_path}.")

        # Return the DataFrame
        return pathways_df

    except Exception as e:
        print(f"An error occurred: {e}")

# Application
gmt_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG_2022_entrez.gmt"
gene_id_dict_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/IDs/entrez2ensembl.txt"
output_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG_2022_entrez_mapped_pathways.tsv"

pathways = process_and_save_pathways(gmt_path, gene_id_dict_path, output_path)
