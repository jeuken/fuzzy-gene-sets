import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def fuzzy_gsea_score(ranked_list_path, pathway_file_path, scores_file_path, plotting_file_path, specific_pathways=None, membership='Default_Membership', plot=False, plot_path=None):
    """
    Calculates enrichment scores for pathways using a ranked list and a dataset of pathways.
    Processes only the specified pathways if specific_pathways is provided, and saves the scores and plotting info to files.
    Optionally generates and saves plots based on the `plot` argument.

    Args:
        ranked_list_path (str): Path to the ranked list file in TSV format.
        pathway_file_path (str): Path to the pathway file in TSV format.
        scores_file_path (str): Path to save the enrichment scores.
        plotting_file_path (str): Path to save the plotting information.
        specific_pathways (list of str or None): List of pathway identifiers to filter by. If None, processes all pathways.
        membership (str): Column name to use for membership values. Options are 'Default_Membership', 'Random_Membership', 'Overlap_Membership'.
        plot (bool): Whether to generate and save plots.
        plot_path (str or None): Directory path where plots should be saved. Required if `plot` is True.
    
    Returns:
        None
    """
    
    # Load the ranked list
    try:
        ranked_df = pd.read_csv(ranked_list_path, sep='\t', header=0)
        if 'Ensembl_ID' not in ranked_df.columns or 'T-statistic' not in ranked_df.columns:
            raise ValueError("Ranked list file must contain 'Ensembl_ID' and 'T-statistic' columns.")
        ranked_df.set_index('Ensembl_ID', inplace=True)
        ranked_list = ranked_df['T-statistic'].to_dict()
    except Exception as e:
        print(f"Error loading ranked list file: {e}")
        return
    
    # Load pathway file
    try:
        pathways_df = pd.read_csv(pathway_file_path, sep='\t')
        if membership not in pathways_df.columns or 'Pathway_Name' not in pathways_df.columns or 'Gene_ID' not in pathways_df.columns:
            raise ValueError(f"Pathway file must contain '{membership}', 'Pathway_Name', and 'Gene_ID' columns.")
    except Exception as e:
        print(f"Error loading pathway file: {e}")
        return
    
    # Filter out rows with NaN values in the membership column
    pathways_df = pathways_df.dropna(subset=[membership])
    
    # Create a dictionary of pathways
    pathways_dict = {}
    for _, row in pathways_df.iterrows():
        pathway_name = row['Pathway_Name']
        description = row['Description']
        gene_id = row['Gene_ID']
        # Use the specified membership column
        membership_value = row[membership]
        
        if pathway_name not in pathways_dict:
            pathways_dict[pathway_name] = {
                'description': description,
                'genes': {}
            }
        
        pathways_dict[pathway_name]['genes'][gene_id] = membership_value
    
    # Determine pathways to process
    if specific_pathways:
        pathway_names = specific_pathways
    else:
        pathway_names = list(pathways_dict.keys())
    
    all_pathways = {name: pathways_dict[name] for name in pathway_names}
    
    # Initialize dictionary to store enrichment scores
    enrichment_scores = {}
    plotting_values_all = {}
    
    # Calculate total number of genes in ranked list
    N = len(ranked_list)
    
    # Convert ranked list keys to a list for indexed access
    ranked_gene_list = list(ranked_list.keys())
    
    # Calculate enrichment scores for each pathway
    for pathway_name, pathway_info in all_pathways.items():
        gene_ids = set(pathway_info['genes'].keys())
        gene_set = gene_ids
        
        # Create correlation dictionary for the pathway genes
        correlation = {gene_id: ranked_list.get(gene_id, 0) for gene_id in gene_set}
        
        # Calculate N_r as the weighted sum of absolute correlations
        N_r = sum(abs(ranked_list.get(gene_id, 0)) * pathway_info['genes'].get(gene_id, 1) for gene_id in gene_set)
        
        N_misses = N - len(gene_ids)
        
        # Initialize variables for enrichment score calculation
        P_hit = 0
        P_miss = 1
        counter = 0
        enrichment_score = 0.0
        plotting_values = []
        
        # Iterate through the ranked list once
        for idx, gene_id in enumerate(ranked_gene_list):
            if gene_id in gene_set:
                # Update P_hit using the membership value
                membership_value = pathway_info['genes'].get(gene_id, 1)  # Use default membership of 1 if not found
                P_hit += (abs(correlation.get(gene_id, 0)* membership_value)/ N_r) 
                counter += 1
            
            # Calculate P_miss
            P_miss = ((idx - counter) + 1) / N_misses
            
            # Update enrichment score if the current score is higher
            if abs(P_hit - P_miss) > abs(enrichment_score):
                enrichment_score = P_hit - P_miss
            
            # Track the enrichment score for plotting
            plotting_values.append(P_hit - P_miss)
        
        # Store the enrichment score and plotting values for the pathway
        enrichment_scores[pathway_name] = {
            'Enrichment_Score': enrichment_score,
            'Description': pathway_info['description']
        }
        plotting_values_all[pathway_name] = plotting_values

    # Convert enrichment scores to DataFrame and include descriptions
    scores_list = [{'Pathway': pathway, 
                    'Enrichment_Score': info['Enrichment_Score'],
                    'Description': info['Description']} 
                   for pathway, info in enrichment_scores.items()]
    
    # Update file paths to include membership
    scores_file_path = scores_file_path.replace('.tsv', f'_{membership}.tsv')
    plotting_file_path = plotting_file_path.replace('.tsv', f'_{membership}.tsv')
    
    scores_df = pd.DataFrame(scores_list)
    scores_df.to_csv(scores_file_path, sep='\t', index=False)
    
    # Save plotting values to file including description
    plotting_values_df = pd.DataFrame([(pathway, idx, value, pathways_dict[pathway]['description']) 
                                        for pathway, values in plotting_values_all.items() 
                                        for idx, value in enumerate(values)],
                                       columns=['Pathway', 'Rank', 'Enrichment_Value', 'Description'])
    plotting_values_df.to_csv(plotting_file_path, sep='\t', index=False)

    # Generate and save plots if required
    if plot and plot_path:
        # Ensure the plot directory exists
        membership_path = os.path.join(plot_path, membership)
        if not os.path.exists(membership_path):
            os.makedirs(membership_path)
        
        # Plot for each pathway
        for pathway, pathway_data in plotting_values_df.groupby('Pathway'):
            plt.figure(figsize=(10, 6))
            plt.plot(pathway_data['Rank'], pathway_data['Enrichment_Value'], label='Enrichment Value')
            plt.xlabel('Rank')
            plt.ylabel('Enrichment Value')
            plt.title(f'Enrichment Plot for {pathway}')
            plt.legend()
            plt.grid(True)
            plt.savefig(os.path.join(membership_path, f'{pathway}_enrichment_plot.png'))
            plt.close()
                     
                    
def shuffle_phenotypes(df):
    shuffled_df = df.copy()
    for col in shuffled_df.columns:
        if col.startswith('healthy') or col.startswith('cancer'):
            shuffled_df[col] = np.random.permutation(shuffled_df[col])
    return shuffled_df

def compute_null_distribution(ranked_list_path, pathway_file_path, n_permutations, membership='Default_Membership'):
    null_distribution = []
    
    ranked_df = pd.read_csv(ranked_list_path, sep='\t', header=0, index_col='Ensembl_ID')
    
    for _ in range(n_permutations):
        shuffled_df = shuffle_phenotypes(ranked_df)
        shuffled_ranked_list_path = ranked_list_path.replace('.txt', '_shuffled.txt')
        shuffled_df.to_csv(shuffled_ranked_list_path, sep='\t', index=False)
        
        fuzzy_gsea_score(shuffled_ranked_list_path, pathway_file_path, 
                         scores_file_path=f'/tmp/enrichment_scores_permutation_{_}.tsv', 
                         plotting_file_path=f'/tmp/plotting_values_permutation_{_}.tsv',
                         membership=membership)
        
        scores_df = pd.read_csv(f'/tmp/enrichment_scores_permutation_{_}.tsv', sep='\t')
        null_distribution.extend(scores_df['Enrichment_Score'].tolist())
    
    return null_distribution

def calculate_p_value(observed_score, null_distribution):
    null_distribution = np.array(null_distribution)
    p_value = np.mean(null_distribution >= observed_score)
    return p_value

def plot_histogram(null_distribution, observed_score, output_path):
    plt.figure(figsize=(10, 6))
    plt.hist(null_distribution, bins=30, edgecolor='k', alpha=0.7, color='grey')
    plt.axvline(x=observed_score, color='r', linestyle='--', label=f'Observed Enrichment Score: {observed_score:.2f}')
    plt.xlabel('Enrichment Score')
    plt.ylabel('Frequency')
    plt.title('Null Distribution of Enrichment Scores')
    plt.legend()
    plt.grid(True)
    plt.savefig(output_path)
    plt.close()

def main():
    # Number of permutations
    n_permutations = 10
    
    # Membership column
    membership = 'Default_Membership'
    
    # Define file paths
    ranked_list_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/gemma/ranked_list.txt"
    pathway_file_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/pathways_with_membership.tsv"
    scores_file_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/GSEA_output/observed_enrichment_scores.tsv"
    plotting_file_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/GSEA_output/observed_plotting_values.tsv"
    p_values_file_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/GSEA_output/enrichment_p_values.tsv"
    plot_output_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/GSEA_output/null_distribution_histogram.png"

    # Calculate observed enrichment scores
    fuzzy_gsea_score(ranked_list_path, pathway_file_path, 
                     scores_file_path=scores_file_path,
                     plotting_file_path=plotting_file_path,
                     membership=membership, 
                     plot=True, 
                     plot_path=os.path.dirname(plot_output_path))
    
    # Compute null distribution
    null_distribution = compute_null_distribution(ranked_list_path, pathway_file_path, n_permutations, membership=membership)
    
    # Read the observed enrichment scores
    observed_scores_df = pd.read_csv(scores_file_path, sep='\t')
    
    # Calculate p-values
    p_values_list = []
    for _, row in observed_scores_df.iterrows():
        pathway = row['Pathway']
        observed_score = row['Enrichment_Score']
        
        p_value = calculate_p_value(observed_score, null_distribution)
        p_values_list.append({
            'Pathway': pathway,
            'Observed_Score': observed_score,
            'P_Value': p_value
        })
        print(f'Pathway: {pathway}, Observed Score: {observed_score:.2f}, p-value: {p_value:.4f}')
    
    # Save p-values to file
    p_values_df = pd.DataFrame(p_values_list)
    p_values_df.to_csv(p_values_file_path, sep='\t', index=False)
    
    # Plot histograms for each observed score
    for _, row in observed_scores_df.iterrows():
        pathway = row['Pathway']
        observed_score = row['Enrichment_Score']
        
        pathway_plot_path = plot_output_path.replace('.png', f'_{pathway}.png')
        plot_histogram(null_distribution, observed_score, pathway_plot_path)


if __name__ == "__main__":
    main()
