import pandas as pd
import matplotlib.pyplot as plt



def fuzzy_enrichment_score(ranked_list_path, pathway_file_path, scores_file_path, plotting_file_path, specific_pathway=None):
    """
    Calculate enrichment scores for pathways using a ranked list and a pathway file.
    Processes only the specified pathway or the first 5 pathways for analysis if not specified,
    and saves the scores and plotting info to files.

    Args:
        ranked_list_path (str): Path to the ranked list file in TSV format.
        pathway_file_path (str): Path to the pathway file in TSV format.
        scores_file_path (str): Path to save the enrichment scores.
        plotting_file_path (str): Path to save the plotting information.
        specific_pathway (str): Pathway identifier to filter by. If None, processes the first 5 pathways.
    
    Returns:
        None
    """
    
    # Load the ranked list
    ranked_df = pd.read_csv(ranked_list_path, sep='\t', header=0, index_col='Ensembl_ID')
    ranked_list = ranked_df['T_stat'].to_dict()
    
    # Load pathway file and create a dictionary of pathways
    pathways_df = pd.read_csv(pathway_file_path, sep='\t')
    pathways_dict = {}
    
    for _, row in pathways_df.iterrows():
        pathway_name = row['Pathway_Name']
        description = row['Description']
        gene_id = row['Gene_ID']
        membership = row['Membership']
        
        if pathway_name not in pathways_dict:
            pathways_dict[pathway_name] = {
                'description': description,
                'genes': {}
            }
        
        pathways_dict[pathway_name]['genes'][gene_id] = membership
    
    # Filter pathways based on the specific pathway or limit to the first 5
    if specific_pathway:
        pathway_names = [specific_pathway]
    else:
        pathway_names = list(pathways_dict.keys())[:5]
    
    first_5_pathways = {name: pathways_dict[name] for name in pathway_names}
    
    # Initialize dictionary to store enrichment scores
    enrichment_scores = {}
    plotting_values_all = {}
    
    # Calculate total number of genes in ranked list
    N = len(ranked_list)
    
    # Convert ranked list keys to a list for indexed access
    ranked_gene_list = list(ranked_list.keys())
    
    # Calculate enrichment scores for each pathway
    for pathway_name, pathway_info in first_5_pathways.items():
        gene_ids = set(pathway_info['genes'].keys())
        gene_set = gene_ids
        
        # Create correlation dictionary for the pathway genes
        correlation = {gene_id: ranked_list.get(gene_id, 0) for gene_id in gene_set}
        
        # Sum of absolute values of T-statistics for genes in the pathway
        N_r = sum(abs(value) for value in correlation.values())
        
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
                # Update P_hit using membership value
                membership_value = pathway_info['genes'].get(gene_id, 1)
                P_hit += abs(correlation.get(gene_id, 0))/ N_r
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
    
    scores_df = pd.DataFrame(scores_list)
    scores_df.to_csv(scores_file_path, sep='\t', index=False)
    
    # Save plotting values to file including description
    plotting_values_df = pd.DataFrame([(pathway, idx, value, pathways_dict[pathway]['description']) 
                                        for pathway, values in plotting_values_all.items() 
                                        for idx, value in enumerate(values)],
                                       columns=['Pathway', 'Rank', 'Enrichment_Value', 'Description'])
    plotting_values_df.to_csv(plotting_file_path, sep='\t', index=False)


def plot_gsea(plotting_file_path, plot_path):
    """
    Plot GSEA enrichment plots using the plotting information file.

    Args:
        plotting_file_path (str): Path to the file containing plotting information.
        plot_path (str): Directory path where plots should be saved.
    
    Returns:
        None
    """
    
    # Load plotting information
    plotting_df = pd.read_csv(plotting_file_path, sep='\t')
    
    # Plot for each pathway
    for pathway, pathway_data in plotting_df.groupby('Pathway'):
        plt.figure(figsize=(10, 6))
        plt.plot(pathway_data['Rank'], pathway_data['Enrichment_Value'], label='Enrichment Score')
        
        # Add a line indicating the final enrichment score
        final_score = pathway_data.loc[pathway_data['Enrichment_Value'].abs().idxmax(), 'Enrichment_Value']
        plt.axhline(y=final_score, color='r', linestyle='--', label=f'Final Enrichment Score: {final_score:.2f}')
        
        # Include description in the title
        description = pathway_data['Description'].iloc[0]
        plt.title(f'GSEA Plot for {pathway}\n{description}')
        plt.xlabel('Rank')
        plt.ylabel('Enrichment Score')
        plt.legend()
        plt.grid(True)
        
        # Save the plot to the specified path
        plt.savefig(f'{plot_path}/{pathway}_gsea_plot.png')
        plt.close()


# Example usage
ranked_list_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/gemma/ranked_list.txt"
pathway_file_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/mapped_pathways.tsv"
scores_file_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/enrichment_scores.tsv"
plotting_file_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/plotting_values.tsv"
specific_pathway = "hsa05211"
plot_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/plots"
# Get the enrichment scores and save to files
fuzzy_enrichment_score(ranked_list_path, pathway_file_path, scores_file_path, plotting_file_path)

# Plot GSEA plots using the plotting information file
plot_gsea(plotting_file_path,plot_path)



# Get the enrichment scores for the specific pathway and save to files
fuzzy_enrichment_score(ranked_list_path, pathway_file_path, scores_file_path, plotting_file_path, specific_pathway)

# Plot GSEA plots using the plotting information file
plot_gsea(plotting_file_path,plot_path)


