#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd

def fuzzy_ora(membership_file, pathway_file, specific_pathway=None, 
              query_membership_type='Membership_Score', pathway_membership_type='Default_Membership', 
              output_file='ora_enrichment_scores.txt'):
    """
    Computes fuzzy intersection sizes between a query set and pathways, and saves the results.

    Args:
        membership_file (str): Path to the query set membership file.
        pathway_file (str): Path to the pathway information file.
        specific_pathway (str, optional): Specific pathway to focus on, or None for all pathways.
        query_membership_type (str, optional): Column name for query set membership.
        pathway_membership_type (str, optional): Column name for pathway membership.
        output_file (str, optional): Path to save the enrichment scores.

    Returns:
        None
    """
    try:
        # Load query set memberships
        query_set = pd.read_csv(membership_file, sep='\t', index_col='Ensembl_ID')[query_membership_type].astype(float).to_dict()

        # Load and preprocess pathway data
        pathways_df = pd.read_csv(pathway_file, sep='\t', dtype={pathway_membership_type: float, 'Gene_ID': str})
        
        # Group pathways and aggregate Gene_IDs and memberships
        pathways_grouped = pathways_df.groupby('Pathway_Name').agg(
            Gene_ID=('Gene_ID', list), 
            Memberships=(pathway_membership_type, list)
        ).reset_index()

        # Filter to a specific pathway if specified
        if specific_pathway:
            pathways_grouped = pathways_grouped[pathways_grouped['Pathway_Name'] == specific_pathway]
        
        # Compute fuzzy intersection sizes
        enrichment_scores = {
            pathway: sum(min(query_set.get(gene, 0), mem) for gene, mem in zip(genes, memberships))
            for pathway, genes, memberships in zip(
                pathways_grouped['Pathway_Name'], pathways_grouped['Gene_ID'], pathways_grouped['Memberships']
            )
        }

        # Save results to a file
        pd.DataFrame(list(enrichment_scores.items()), columns=['Pathway_Name', 'Fuzzy_Intersection_Size']).to_csv(
            output_file, sep='\t', index=False
        )

        print(f"Enrichment scores saved to {output_file}.")
        
    except Exception as e:
        print(f"Error: {e}")

# Example usage
membership_file = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/gemma/membership_scores.txt"
pathway_file = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/mapped_pathways.tsv"
specific_pathway = None  # Set to a specific pathway name or leave as None for all pathways
query_membership_type = "Membership_Score"
pathway_membership_type = "Overlap_Membership"
output_file = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/ORA_output/ora_enrichment_scores.txt"

# Call the function with all required arguments
fuzzy_ora(membership_file, pathway_file, specific_pathway, query_membership_type, pathway_membership_type, output_file)
