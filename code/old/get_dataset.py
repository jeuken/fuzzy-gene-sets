import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
import os

def shuffle_phenotype_labels(expression_df):
    """
    Shuffle the phenotype labels in the expression DataFrame.
    
    Args:
        expression_df (pd.DataFrame): DataFrame with phenotype labels to shuffle.
    
    Returns:
        pd.DataFrame: DataFrame with shuffled phenotype labels.
    """
    shuffled_df = expression_df.copy()
    columns = shuffled_df.columns
    shuffled_df.columns = np.random.permutation(columns)
    return shuffled_df

def create_ranked_list(expression_df):
    """
    Create a ranked list of genes based on T-statistics from the expression DataFrame.
    
    Args:
        expression_df (pd.DataFrame): DataFrame containing expression data with phenotype labels.
    
    Returns:
        dict: Dictionary with gene IDs as keys and T-statistics as values.
    """
    healthy_cols = [col for col in expression_df.columns if col.startswith('healthy')]
    cancer_cols = [col for col in expression_df.columns if col.startswith('cancer')]
    
    healthy_data = expression_df[healthy_cols].values
    cancer_data = expression_df[cancer_cols].values
    
    t_stats, _ = ttest_ind(healthy_data, cancer_data, axis=1, equal_var=False, nan_policy='omit')
    
    ranked_list = dict(zip(expression_df['Ensembl_ID'], t_stats))
    ranked_list = {k: v for k, v in sorted(ranked_list.items(), key=lambda item: item[1], reverse=True)}
    return ranked_list

def compute_fuzzy_enrichment_score(ranked_list, pathway_info):
    """
    Compute the fuzzy enrichment score for a pathway based on the ranked list of genes.
    
    Args:
        ranked_list (dict): Dictionary with gene IDs as keys and T-statistics as values.
        pathway_info (dict): Dictionary with pathway gene information.
    
    Returns:
        float: Fuzzy enrichment score.
        list: List of enrichment values for plotting.
    """
    gene_ids = set(pathway_info['genes'].keys())
    N_r = sum(abs(ranked_list.get(gene_id, 0)) * pathway_info['genes'].get(gene_id, 1) for gene_id in gene_ids)
    N = len(ranked_list)
    N_misses = N - len(gene_ids)
    
    P_hit = 0
    P_miss = 1
    counter = 0
    enrichment_score = 0.0
    plotting_values = []
    
    for idx, gene_id in enumerate(ranked_list):
        if gene_id in gene_ids:
            membership_value = pathway_info['genes'].get(gene_id, 1)
            P_hit += abs(ranked_list.get(gene_id, 0)) * membership_value / N_r
            counter += 1
        P_miss = ((idx - counter) + 1) / N_misses
        enrichment_score = max(enrichment_score, P_hit - P_miss)
        plotting_values.append(P_hit - P_miss)
    
    return enrichment_score, plotting_values

def make_enrichment_score_plot(plotting_values, pathway_name, plot_path, membership, description):
    """
    Create and save the enrichment score plot for a pathway.
    
    Args:
        plotting_values (list): List of enrichment values for plotting.
        pathway_name (str): Name of the pathway.
        plot_path (str): Directory path to save the plot.
        membership (str): Membership type.
        description (str): Description of the pathway.
    """
    if not os.path.exists(plot_path):
        os.makedirs(plot_path)
    
    plt.figure(figsize=(10, 6))
    plt.plot(plotting_values)
    plt.title(f'GSEA Plot for {pathway_name}\nMembership: {membership}\n{description}')
    plt.xlabel('Rank')
    plt.ylabel('Enrichment Score')
    plt.grid(True)
    plt.savefig(os.path.join(plot_path, f'{pathway_name}_{membership}_gsea_plot.png'))
    plt.close()

def create_null_distribution(expression_df, ranked_list, pathway_info, num_permutations=1000):
    """
    Create a null distribution of enrichment scores by permuting phenotype labels.
    
    Args:
        expression_df (pd.DataFrame): DataFrame with expression data.
        ranked_list (dict): Dictionary with gene IDs as keys and T-statistics as values.
        pathway_info (dict): Dictionary with pathway gene information.
        num_permutations (int): Number of permutations to perform.
    
    Returns:
        list: Null distribution of enrichment scores.
    """
    null_distribution = []
    
    for _ in range(num_permutations):
        shuffled_df = shuffle_phenotype_labels(expression_df)
        permuted_ranked_list = create_ranked_list(shuffled_df)
        enrichment_score, _ = compute_fuzzy_enrichment_score(permuted_ranked_list, pathway_info)
        null_distribution.append(enrichment_score)
    
    return null_distribution

def main(expression_path, pathway_file_path, scores_file_path, plot_path, num_permutations=1000, membership='Default_Membership'):
    """
    Main function to compute observed enrichment scores, compare to null distribution, 
    and save plots and results.
    
    Args:
        expression_path (str): Path to the expression data file.
        pathway_file_path (str): Path to the pathway file.
        scores_file_path (str): Path to save the enrichment scores.
        plot_path (str): Directory path to save the plots.
        num_permutations (int): Number of permutations for null distribution.
        membership (str): Pathway membership type to use.
    """
    # Load expression data
    expression_df = pd.read_csv(expression_path, sep='\t')
    
    # Load pathway data
    pathways_df = pd.read_csv(pathway_file_path, sep='\t')
    pathways_dict = {}
    
    for _, row in pathways_df.iterrows():
        pathway_name = row['Pathway_Name']
        description = row['Description']
        gene_id = row['Gene_ID']
        membership_value = row.get(membership, 1)
        
        if pathway_name not in pathways_dict:
            pathways_dict[pathway_name] = {
                'description': description,
                'genes': {}
            }
        
        pathways_dict[pathway_name]['genes'][gene_id] = membership_value
    
    # Compute observed enrichment scores and save results
    enrichment_scores = {}
    for pathway_name, pathway_info in pathways_dict.items():
        ranked_list = create_ranked_list(expression_df)
        enrichment_score, plotting_values = compute_fuzzy_enrichment_score(ranked_list, pathway_info)
        
        enrichment_scores[pathway_name] = enrichment_score
        make_enrichment_score_plot(plotting_values, pathway_name, plot_path, membership, pathway_info['description'])
    
    # Create and save null distribution
    null_distribution = {pathway: create_null_distribution(expression_df, create_ranked_list(expression_df), pathway_info, num_permutations) for pathway, pathway_info in pathways_dict.items()}
    
    # Compute p-values and save results
    p_values = {pathway: np.mean(np.array(null_distribution[pathway]) >= score) for pathway, score in enrichment_scores.items()}
    
    # Save enrichment scores with p-values
    scores_df = pd.DataFrame({
        'Pathway': list(enrichment_scores.keys()),
        'Enrichment_Score': list(enrichment_scores.values()),
        'p-value': [p_values[pathway] for pathway in enrichment_scores.keys()]
    })
    scores_df.to_csv(scores_file_path, sep='\t', index=False)
    
    # Save null distribution histograms
    for pathway, scores in null_distribution.items():
        plt.figure(figsize=(10, 6))
        plt.hist(scores, bins=30, edgecolor='k')
        plt.title(f'Null Distribution for {pathway}')
        plt.xlabel('Enrichment Score')
        plt.ylabel('Frequency')
        plt.savefig(os.path.join(plot_path, f'{pathway}_null_distribution.png'))
        plt.close()

    print(f"Enrichment scores and null distributions saved to {scores_file_path} and {plot_path}.")

# Example usage
expression_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/gemma/expression_dataframe.txt"
pathway_file_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/mapped_pathways.tsv"
scores_file_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/GSEA_output/enrichment_scores.tsv"
plot_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/GSEA_output/plots"
membership = 'Overlap_Membership'

main(expression_path, pathway_file_path, scores_file_path, plot_path, num_permutations=10, membership=membership)
