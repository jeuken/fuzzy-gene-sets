import pandas as pd
import matplotlib.pyplot as plt

# Paths to your files
ensembl_file_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/MalaCards_Ensembl.csv"
pathway_file_path = "../../data/pathways/KEGG/KEGG_2022_membership.tsv"  # Path to the pathway file
fuzzy_enrichment_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/ORA_output/T2D_E-HCAD-31/ora_results_T2D_E-HCAD-31_Crisp_Membership_Overlap_Membership_linear.csv"  # Fuzzy enrichment results file
standard_enrichment_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/ORA_output/T2D_E-HCAD-31/standard_ora_results_T2D_E-HCAD-31_Crisp_Membership_Crisp_Membership.csv"  # Standard enrichment results file

# Import the Ensembl conversion file
ensembl_df = pd.read_csv(ensembl_file_path)

# Import the pathway file (specifying that it's tab-separated)
pathway_df = pd.read_csv(pathway_file_path, sep='\t')  # Use sep='\t' for TSV files

# Merge the two DataFrames on 'Ensembl_ID'
merged_df = pd.merge(pathway_df, ensembl_df[['Ensembl_ID', 'Score', 'Symbol']], on='Ensembl_ID', how='left')

# Replace NaN scores with 0 (without inplace)
merged_df['Score'] = merged_df['Score'].fillna(0)

# Calculate the average score for each pathway, including Description
average_scores = merged_df.groupby(['Pathway_Name', 'Description'])['Score'].mean().reset_index()

# Function to collect genes information
def get_genes_info(group):
    non_zero_genes = group[group['Score'] > 0]
    ensembl_ids = non_zero_genes['Ensembl_ID'].tolist()
    symbols = non_zero_genes['Symbol'].tolist()
    gene_count = non_zero_genes['Ensembl_ID'].nunique()
    return pd.Series({
        'Genes': ', '.join(ensembl_ids),
        'Symbols': ', '.join(symbols),
        'Gene_Count': gene_count
    })

# Group by Pathway_Name and Description to get the genes information
genes_info = merged_df.groupby(['Pathway_Name', 'Description'], as_index=False, group_keys=False).apply(get_genes_info)

# Merge the genes info back into the average_scores DataFrame
average_scores = pd.merge(average_scores, genes_info, on=['Pathway_Name', 'Description'])

# Filter for pathways with at least 4 genes
average_scores = average_scores[average_scores['Gene_Count'] >= 4]

# Function to calculate total weighted score for an enrichment file
def calculate_weighted_scores(enrichment_results_path, label):
    # Load enrichment results
    enrichment_results = pd.read_csv(enrichment_results_path)
    
    # Merge enrichment results with average_scores
    result_df = pd.merge(average_scores, enrichment_results[['Pathway_Name', 'Rank']], on='Pathway_Name', how='inner')
    
    # Calculate the weight
    total_pathways = len(enrichment_results)
    result_df['Weight'] = 1 - (result_df['Rank'] / total_pathways)
    
    # Calculate the weighted score
    result_df['Weighted_Score'] = result_df['Weight'] * result_df['Score']
    
    # Sum the total weighted score
    total_weighted_score = result_df['Weighted_Score'].sum()
    
    # Rename columns for clarity
    result_df.rename(columns={
        'Rank': f'Rank_{label}',
        'Weight': f'Weight_{label}',
        'Weighted_Score': f'Weighted_Score_{label}'
    }, inplace=True)
    
    return result_df, total_weighted_score

# Calculate weighted scores for both fuzzy and standard enrichment results
fuzzy_results_df, fuzzy_total_weighted_score = calculate_weighted_scores(fuzzy_enrichment_path, 'Fuzzy')
standard_results_df, standard_total_weighted_score = calculate_weighted_scores(standard_enrichment_path, 'Standard')

# Combine fuzzy and standard results into one DataFrame
combined_results_df = pd.merge(
    fuzzy_results_df,
    standard_results_df[['Pathway_Name', 'Rank_Standard', 'Weight_Standard','Weighted_Score_Standard']],
    on='Pathway_Name',
    how='outer'
)

# Display the total weighted scores for both analyses
print(f"Total Weighted Score (Fuzzy): {fuzzy_total_weighted_score}")
print(f"Total Weighted Score (Standard): {standard_total_weighted_score}")

# Optional: Display the combined results
print(combined_results_df[['Pathway_Name', 'Description', 'Score', 'Rank_Fuzzy', 'Weight_Fuzzy', 'Weighted_Score_Fuzzy', 'Rank_Standard', 'Weight_Standard', 'Weighted_Score_Standard']])

# Scatter plots
plt.figure(figsize=(12, 5))

# Scatter plot for Score vs Rank_Fuzzy
plt.subplot(1, 2, 1)
plt.scatter(combined_results_df['Rank_Fuzzy'], combined_results_df['Score'], color='blue', alpha=0.7)
plt.title('Score vs Rank (Fuzzy)')
plt.xlabel('Rank (Fuzzy)')
plt.ylabel('Score')
plt.grid(True)

# Scatter plot for Score vs Rank_Standard
plt.subplot(1, 2, 2)
plt.scatter(combined_results_df['Rank_Standard'], combined_results_df['Score'], color='green', alpha=0.7)
plt.title('Score vs Rank (Standard)')
plt.xlabel('Rank (Standard)')
plt.ylabel('Score')
plt.grid(True)

# Adjust layout and show the plots
plt.tight_layout()
plt.show()
