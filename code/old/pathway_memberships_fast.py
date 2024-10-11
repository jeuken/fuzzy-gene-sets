#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import os
from tqdm import tqdm

def process_and_expand_pathways(gmt_path, ensembl_dict_path, string_scores_path, output_file):
    """
    Process and expand pathways using GMT file, Ensembl dictionary, and STRING scores.
    """
    gene_dict = pd.read_csv(ensembl_dict_path, sep='\t', header=1, 
                             names=['Ensembl_ID', 'Gene_Name'], 
                             dtype={'Gene_Name': str}
                            ).dropna().drop_duplicates(subset=['Gene_Name']).set_index('Gene_Name')

    pathway_rows, conversion_errors = [], {}

    with open(gmt_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                pathway_name, description, *gene_names = parts
                gene_names = [name.strip() for name in gene_names]
                gene_names = list(set(gene_names))  # Remove duplicate gene names for the pathway
                total_genes = len(gene_names)
                valid_gene_names = [name for name in gene_names if name in gene_dict.index]
                invalid_gene_count = total_genes - len(valid_gene_names)
                
                if valid_gene_names:
                    ensembl_ids = gene_dict.loc[valid_gene_names, 'Ensembl_ID'].values
                    ensembl_ids = list(set(ensembl_ids))  # Remove duplicate Ensembl IDs for the pathway
                    pathway_rows.extend([[pathway_name, description, ensembl_id] for ensembl_id in ensembl_ids])

                conversion_errors[pathway_name] = {
                    'Number_Not_Converted': invalid_gene_count,
                    'Total_Genes': total_genes,
                    'Percentage_Not_Converted': (invalid_gene_count / total_genes) * 100 if total_genes > 0 else 0
                }

    pathways_df = pd.DataFrame(pathway_rows, columns=['Pathway_Name', 'Description', 'Ensembl_ID'])
    pathways_df['Crisp_Membership'] = 1
    pathways_df['Expansion_Membership'] = 1
    
    pathway_conversion_errors = pd.DataFrame.from_dict(conversion_errors, orient='index').reset_index()
    pathway_conversion_errors.columns = ['Pathway_Name', 'Number_Not_Converted', 'Total_Genes', 'Percentage_Not_Converted']
    pathway_conversion_errors['Percentage_Not_Converted'] = pathway_conversion_errors['Percentage_Not_Converted'].round(2)

    gene_counts = pathways_df['Ensembl_ID'].value_counts()
    max_freq = gene_counts.max()
    min_freq = gene_counts.min()
    cumulative_distribution = gene_counts.rank(pct=True, method='max',ascending = False)
    print("min",cumulative_distribution.min())
    print("may",cumulative_distribution.max())
    pathways_df['Overlap_Membership'] = cumulative_distribution[pathways_df['Ensembl_ID']].values
    pathways_df['Overlap_Membership_linear'] = (max_freq - pathways_df['Ensembl_ID'].map(gene_counts)) / (max_freq - min_freq)

    # Load STRING scores
    string_scores = pd.read_csv(string_scores_path, sep="\t")
    
    def process_pathway(pathway_name, pathway_df):
        pathway_genes = set(pathway_df['Ensembl_ID'])
        gene1_scores = string_scores[string_scores['Gene2'].isin(pathway_genes)]
        gene1_scores = gene1_scores[gene1_scores['string_score'] >= 0.4]
        
        if gene1_scores.empty:
            return pd.DataFrame()  # Return an empty DataFrame if no genes are left

        gene1_scores = gene1_scores.groupby('Gene1')['string_score'].mean().reset_index()
        percentiles = gene1_scores['string_score'].rank(pct=True).values
        gene_percentile_map = dict(zip(gene1_scores['Gene1'], percentiles))
    
        return pd.DataFrame([{
            'Pathway_Name': pathway_name,
            'Ensembl_ID': gene,
            'Expansion_Membership': gene_percentile_map.get(gene, 0),
            'Description': pathway_df['Description'].iloc[0],
            'Crisp_Membership': np.nan,
            'Overlap_Membership': np.nan,
            
        } for gene in gene1_scores['Gene1'] if gene not in pathway_genes])

    # Adding progress bar for processing pathways
    new_genes_df = pd.concat(
        [process_pathway(name, group) for name, group in tqdm(pathways_df.groupby('Pathway_Name'), desc="Processing Pathways")],
        ignore_index=True
    )

    final_df = pd.concat([pathways_df, new_genes_df], ignore_index=True)
    final_df.to_csv(output_file, sep="\t", index=False)
    print(f"Filtered and expanded DataFrame saved to {output_file}")

    percentage_table_path = os.path.join(os.path.dirname(output_file), 'percentage_not_matched.tsv')
    pathway_conversion_errors.to_csv(percentage_table_path, sep='\t', index=False)
    print(f"Percentage of genes not matched saved to {percentage_table_path}")

    return final_df


# Example usage
gmt_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/KEGG_2022_entrez.gmt"
ensembl_dict_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/IDs/entrez2ensembl.txt"
string_scores_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/STRING/STRING_GGI.txt"
output_file = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/KEGG_2022_entrez_with_membership_new.tsv"

final_df = process_and_expand_pathways(gmt_path, ensembl_dict_path, string_scores_path, output_file)
