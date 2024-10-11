import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from scipy.stats import ttest_ind

# Shuffle phenotype labels
def shuffle_phenotype_labels(df):
    shuffled_df = df.copy()
    if shuffled_df.index.name == 'Ensembl_ID':
        shuffled_df.reset_index(inplace=True)

    phenotype_cols = [col for col in shuffled_df.columns if col.startswith(('control', 'disease'))]
    shuffled_col_names = np.random.permutation(phenotype_cols)
    shuffled_df = shuffled_df.rename(columns=dict(zip(phenotype_cols, shuffled_col_names)))

    return shuffled_df

# Calculate statistics
def calculate_statistics(df):
    if 'Ensembl_ID' in df.columns:
        df.set_index('Ensembl_ID', inplace=True)

    control_cols = [col for col in df.columns if col.startswith('control')]
    disease_cols = [col for col in df.columns if col.startswith('disease')]

    df[control_cols + disease_cols] = df[control_cols + disease_cols].apply(pd.to_numeric, errors='coerce')
    df.dropna(subset=control_cols + disease_cols, inplace=True)

    control_data, disease_data = df[control_cols].values, df[disease_cols].values
    _, p_vals = ttest_ind(control_data, disease_data, axis=1, equal_var=False, nan_policy='omit')

    adjusted_p_vals = multipletests(p_vals, method='fdr_bh')[1]
    return pd.DataFrame({'Ensembl_ID': df.index, 'adj_p_value': adjusted_p_vals}).reset_index(drop=True)

# Compute query membership based on adjusted p-values
def compute_query_membership(adj_p_value, membership_type='fuzzy'):
    if membership_type == 'crisp':
        return 0 if adj_p_value > 0.05 else 1
    else:
        return 1 - adj_p_value

# Compute fuzzy intersection size between query set and a pathway
def compute_fuzzy_intersection(query_set, pathway_genes, pathway_memberships):
    intersection_size = 0
    for gene, membership in zip(pathway_genes, pathway_memberships):
        query_value = query_set.get(gene, 0)
        min_value = min(query_value, membership)
        intersection_size += min_value
        print(f"Gene: {gene}, Pathway Membership: {membership}, Query Membership: {query_value}, min(Query, Pathway): {min_value}, Accumulated Intersection: {intersection_size}")
    return intersection_size

# Permute and calculate null distribution
def permute_and_calculate_null_distribution(expression_df, pathways_grouped, query_membership, n_permutations=1000):
    null_distributions = {pathway: [] for pathway in pathways_grouped['Pathway_Name']}

    for i in range(n_permutations):
        shuffled_df = shuffle_phenotype_labels(expression_df)
        ranked_df = calculate_statistics(shuffled_df)
        query_set = {gene: compute_query_membership(p_val, query_membership) for gene, p_val in ranked_df.set_index('Ensembl_ID')['adj_p_value'].items()}

        for pathway_name, genes, memberships in zip(
                pathways_grouped['Pathway_Name'], pathways_grouped['Ensembl_ID'], pathways_grouped['Memberships']):
            score = compute_fuzzy_intersection(query_set, genes, memberships)
            null_distributions[pathway_name].append(score)

        if (i + 1) % 100 == 0:
            print(f'Completed {i + 1} permutations.')

    return null_distributions

# Calculate one-sided p-value for observed scores
def calculate_p_value(observed_score, null_distribution):
    return sum(1 for score in null_distribution if score >= observed_score) / len(null_distribution)

# Save null distribution plots
def save_null_distribution_plots(null_distributions, observed_scores, pathways_grouped, plot_path, membership):
    plot_dir = os.path.join(plot_path, membership, 'null_distributions')
    os.makedirs(plot_dir, exist_ok=True)

    for pathway in null_distributions:
        plt.figure(figsize=(10, 6))
        plt.hist(null_distributions[pathway], bins=30, alpha=0.75, color='blue', label='Null Distribution')
        plt.axvline(observed_scores[pathway], color='red', linestyle='--', linewidth=2, label=f'Observed: {observed_scores[pathway]:.2f}')
        plt.title(f'Null Distribution for {pathway}')
        plt.xlabel('Enrichment Score')
        plt.ylabel('Frequency')
        plt.legend()
        plt.grid(True)
        plt.savefig(os.path.join(plot_dir, f'{pathway}_{membership}_null_distribution.png'))
        plt.close()

# Main function for gene set enrichment analysis
def main(expression_path, pathway_file_path, output_path, plot_path=None, membership='Overlap_Membership', query_membership='fuzzy', n_permutations=1000):
    try:
        print("Loading expression data...")
        expression_df = pd.read_csv(expression_path)
        print(f"Expression data loaded with {expression_df.shape[0]} genes and {expression_df.shape[1]} columns.")

        print("Calculating statistics...")
        ranked_df = calculate_statistics(expression_df)
        print(f"Statistics calculated. Adjusted p-values (first 10):\n{ranked_df.head(10)}")

        print("Computing query memberships...")
        query_set = {gene: compute_query_membership(p_val, query_membership) for gene, p_val in ranked_df.set_index('Ensembl_ID')['adj_p_value'].items()}
        print(f"Query Set Memberships (first 10): {list(query_set.items())[:10]}")

        print("Loading pathways...")
        pathways_df = pd.read_csv(pathway_file_path, sep='\t').dropna(subset=[membership])

        # Group pathways by 'Pathway_Name' and list their associated genes and memberships
        pathways_grouped = pathways_df.groupby('Pathway_Name').agg(
            {'Ensembl_ID': list, membership: list}
        ).reset_index()
        pathways_grouped.rename(columns={membership: 'Memberships'}, inplace=True)

        print("Computing observed fuzzy intersection sizes...")
        observed_scores = {
            pathway: compute_fuzzy_intersection(query_set, genes, memberships)
            for pathway, genes, memberships in zip(
                pathways_grouped['Pathway_Name'], pathways_grouped['Ensembl_ID'], pathways_grouped['Memberships']
            )
        }

        print(f"Calculating null distributions with {n_permutations} permutations...")
        null_distributions = permute_and_calculate_null_distribution(expression_df, pathways_grouped, query_membership, n_permutations)
        print("Null distributions calculated.")

        print("Saving results...")
        results = []
        for pathway in observed_scores:
            p_value = calculate_p_value(observed_scores[pathway], null_distributions[pathway])
            results.append({'Pathway_Name': pathway, 'Observed_Score': observed_scores[pathway], 'p-value': p_value})

        results_df = pd.DataFrame(results).sort_values(by='p-value')
        results_file_path = os.path.join(output_path, f'{membership}_{query_membership}_results.tsv')
        results_df.to_csv(results_file_path, sep='\t', index=False)
        print(f"Results saved to {results_file_path}.")

        if plot_path:
            print(f"Saving null distribution plots to {plot_path}...")
            save_null_distribution_plots(null_distributions, observed_scores, pathways_grouped, plot_path, query_membership)
            print(f"Plots saved to {plot_path}.")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    expression_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/gemma/Alzheimers_GSE95587_formatted.csv"
    pathway_file_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/KEGG_2022_entrez_with_membership.tsv"
    output_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/ORA_output/Alzheimers_GSE95587_sp"
    plot_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/ORA_output/Alzheimers_GSE95587_sp/null_distributions"
    
    main(expression_path, pathway_file_path, output_path, plot_path=plot_path, membership='Overlap_Membership', query_membership='crisp', n_permutations=100)
