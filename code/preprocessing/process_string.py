#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import os

def process_string(protein2gene_path, string_ppi_path, string_ggi_path):
    """
    Processes the STRING PPI data by mapping protein IDs to gene IDs using a provided gene-protein mapping file.
    The processed data is saved as a gene-gene interaction (GGI) file.

    Args:
        protein2gene_path (str): Path to the protein-to-gene mapping file.
        string_ppi_path (str): Path to the STRING PPI file.
        string_ggi_path (str): Path to save the processed STRING GGI file.

    Returns:
        pd.DataFrame: DataFrame containing the processed gene-gene interaction information.
    """
    try:
        # Load and preprocess the gene-protein mapping
        gene_dict = (
            pd.read_csv(protein2gene_path, sep='\t', dtype=str)
            .assign(
                Gene_ID=lambda df: df['Gene stable ID'].str.strip(),
                Protein_ID=lambda df: df['Protein stable ID'].str.strip()
            )
            .dropna(subset=['Protein_ID'])
            .loc[lambda df: df['Protein_ID'] != '']
            .drop_duplicates(subset=['Protein_ID'])
            .set_index('Protein_ID')['Gene_ID']
            .to_dict()
        )

        print(f"Loaded gene dictionary with {len(gene_dict)} entries.")

        # Check if the STRING PPI file exists
        if not os.path.exists(string_ppi_path):
            raise FileNotFoundError(f"STRING PPI file not found: {string_ppi_path}")

        # Read the STRING PPI file with flexible whitespace handling
        ppi_data = pd.read_csv(string_ppi_path, sep='\s+', dtype=str)

        # Verify the header format
        expected_columns = ['protein1', 'protein2', 'combined_score']
        if list(ppi_data.columns) != expected_columns:
            raise ValueError(f"Unexpected header format in STRING PPI file: {list(ppi_data.columns)}")

        # Clean and map proteins to genes
        ppi_data = (
            ppi_data.assign(
                protein1=lambda df: df['protein1'].str.replace('9606.', '').str.strip(),
                protein2=lambda df: df['protein2'].str.replace('9606.', '').str.strip()
            )
            .assign(
                Gene1=lambda df: df['protein1'].map(gene_dict),
                Gene2=lambda df: df['protein2'].map(gene_dict)
            )
            .dropna(subset=['Gene1', 'Gene2'])  # Drop rows with unmapped genes
            .rename(columns={'combined_score': 'string_score'})[['Gene1', 'Gene2', 'string_score']]  # Select relevant columns
        )

        print(f"Processed STRING data with {ppi_data.shape[0]} interactions after gene mapping.")

        # Save the processed data to a TSV file
        ppi_data.to_csv(string_ggi_path, index=False, sep='\t')
        print(f"Processed STRING GGI data saved to {string_ggi_path}.")

        return ppi_data

    except Exception as e:
        print(f"An error occurred: {e}")

# Example usage
protein2gene_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/IDs/protein2gene.txt"
string_ppi_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/STRING/STRING_PPI.txt"
string_ggi_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/STRING/STRING_GGI.txt"

process_string(protein2gene_path, string_ppi_path, string_ggi_path)

