import pandas as pd
import os
from tqdm import tqdm
import matplotlib.pyplot as plt

def load_dataframe(csv_path):
    return pd.read_csv(csv_path, sep=",")

def add_crisp_membership(pathways_df):
    pathways_df['Crisp_Membership'] = 1
    return pathways_df

def calculate_overlap_membership(pathways_df):
    gene_counts = pathways_df['Ensembl_ID'].value_counts()
    
    # Add Frequency column to the pathways DataFrame
    pathways_df['Frequency'] = pathways_df['Ensembl_ID'].map(gene_counts)
    
    max_freq = gene_counts.max()
    min_freq = gene_counts.min()
    total_pathways = pathways_df['Pathway_Name'].nunique()
    
    # Calculate cumulative distributions
    cumulative_distribution_max = gene_counts.rank(pct=True, method='max', ascending=False)
    #cumulative_distribution_min = gene_counts.rank(pct=True, method='min', ascending=False)
    #cumulative_distribution_av = gene_counts.rank(pct=True, method='average', ascending=False)
    # Calculate cumulative distributions
    #cumulative_distribution_max_a = gene_counts.rank(pct=True, method='max', ascending=True)
    #cumulative_distribution_min_a = gene_counts.rank(pct=True, method='min', ascending=True)
    #cumulative_distribution_av_a = gene_counts.rank(pct=True, method='average', ascending=True)

    pathways_df['Overlap_Membership_max'] = cumulative_distribution_max[pathways_df['Ensembl_ID']].values
    #pathways_df['Overlap_Membership_average'] = cumulative_distribution_av[pathways_df['Ensembl_ID']].values
    #pathways_df['Overlap_Membership_min'] = cumulative_distribution_min[pathways_df['Ensembl_ID']].values
    #pathways_df['Overlap_Membership_max_a'] = 1-cumulative_distribution_max_a[pathways_df['Ensembl_ID']].values
    #pathways_df['Overlap_Membership_average_a'] = 1-cumulative_distribution_av_a[pathways_df['Ensembl_ID']].values
    #pathways_df['Overlap_Membership_min_a'] = 1-cumulative_distribution_min_a[pathways_df['Ensembl_ID']].values
    pathways_df['Overlap_Membership_linear'] = (max_freq  - pathways_df['Ensembl_ID'].map(gene_counts) + 1) / (max_freq - min_freq + 1)
    pathways_df['Overlap_Membership'] = (total_pathways - pathways_df['Ensembl_ID'].map(gene_counts)) / (total_pathways - 1)
    
    return pathways_df

def calculate_expansion_membership(pathways_df, string_scores_path):
    string_scores = pd.read_csv(string_scores_path, sep="\t")

    def expansion_membership(pathway_name, pathway_df):
        pathway_genes = set(pathway_df['Ensembl_ID'])
        gene1_scores = string_scores[string_scores['Gene2'].isin(pathway_genes)]
        gene1_scores = gene1_scores[gene1_scores['string_score'] >= 0.4]

        if gene1_scores.empty:
            return pd.DataFrame()

        gene1_scores = gene1_scores.groupby('Gene1')['string_score'].mean().reset_index()
        percentiles = gene1_scores['string_score'].rank(pct=True).values
        gene_percentile_map = dict(zip(gene1_scores['Gene1'], percentiles))

        return pd.DataFrame([{
            'Pathway_Name': pathway_name,
            'Ensembl_ID': gene,
            'Expansion_Membership': gene_percentile_map.get(gene, 0),
            'Description': pathway_df['Description'].iloc[0],
        } for gene in gene1_scores['Gene1'] if gene not in pathway_genes])

    new_genes_df = pd.concat(
        [expansion_membership(name, group) for name, group in tqdm(pathways_df.groupby('Pathway_Name'), desc="Processing Pathways")],
        ignore_index=True
    )
    return new_genes_df


def plot_membership_histograms(pathways_df, output_path):
    # Create a plots directory if it doesn't exist
    plots_dir = os.path.join(output_path, 'plots')
    os.makedirs(plots_dir, exist_ok=True)

    # Get all membership columns once (outside the loop)
    membership_cols = [col for col in pathways_df.columns if 'Membership' in col]

    # Create separate directories for each membership column once
    for column in membership_cols:
        membership_dir = os.path.join(plots_dir, column)
        os.makedirs(membership_dir, exist_ok=True)

    # Group by pathway to create separate plots
    for pathway_name, group in pathways_df.groupby('Pathway_Name'):
        for column in membership_cols:
            plt.figure(figsize=(10, 6))
            group[column].hist(bins=20, color='skyblue', edgecolor='black', alpha=0.7)
            plt.title(f'Histogram of {column} Values for {pathway_name}', fontsize=14, weight='bold')
            plt.xlabel(column)
            plt.ylabel('Frequency')
            plt.xlim(0, 1)  # Set x-axis limits from 0 to 1
            plt.grid(True, linestyle='--', alpha=0.6)
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, column, f'histogram_{pathway_name}.png'))
            plt.close()

    # Deduplicate genes for overall histograms (only consider each gene once)
    deduplicated_df = pathways_df.groupby('Ensembl_ID').agg({
        **{col: 'max' for col in membership_cols},  # Use max membership for each gene across pathways
        'Frequency': 'first'  # Keep the frequency as is (it should be the same for each gene)
    }).reset_index()

    # Create a single histogram across all pathways for each membership
    for column in membership_cols:
        plt.figure(figsize=(10, 6))
        deduplicated_df[column].hist(bins=30, color='salmon', edgecolor='black', alpha=0.7)
        plt.title(f'Histogram of {column} Values Across All Pathways (Deduplicated)', fontsize=14, weight='bold')
        plt.xlabel(column)
        plt.ylabel('Frequency')
        plt.xlim(0, 1)  # Set x-axis limits from 0 to 1
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, f'histogram_all_{column}_deduplicated.png'))
        plt.close()

    # Plot frequency of gene occurrences across all pathways
    plt.figure(figsize=(10, 6))
    deduplicated_df['Frequency'].hist(bins=30, color='lightgreen', edgecolor='black', alpha=0.7)
    plt.title('Histogram of Gene Occurrences Across All Pathways (Deduplicated)', fontsize=14, weight='bold')
    plt.xlabel('Frequency (Number of Pathways a Gene Appears In)', fontsize=12)
    plt.ylabel('Count of Genes', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'histogram_gene_frequency_deduplicated.png'))
    plt.close()



def plot_membership_scatterplots(pathways_df, output_path):
    # Create a scatterplots directory if it doesn't exist
    scatterplots_dir = os.path.join(output_path, 'scatterplots')
    os.makedirs(scatterplots_dir, exist_ok=True)

    # Get all membership columns once (outside the loop)
    membership_cols = [col for col in pathways_df.columns if 'Membership' in col]

    # Generate scatterplots for each membership column
    for column in membership_cols:
        plt.figure(figsize=(10, 6))
        plt.scatter(pathways_df['Frequency'], pathways_df[column], alpha=0.5, color='blue', edgecolor='k')
        plt.title(f'Scatterplot of {column} vs Gene Frequency', fontsize=14, weight='bold')
        plt.xlabel('Frequency (Number of Pathways a Gene Appears In)', fontsize=12)
        plt.ylabel(f'{column}', fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.tight_layout()
        plt.savefig(os.path.join(scatterplots_dir, f'scatter_{column}_vs_frequency.png'))
        plt.close()

    print(f"Scatterplots saved to {scatterplots_dir}")

def process_and_expand_pathways(dataframe_path, string_scores_path, output_file, memberships, plot):
    pathways_df = load_dataframe(dataframe_path)
    print(pathways_df.head())

    if 'crisp' in memberships:
        pathways_df = add_crisp_membership(pathways_df)
    if 'overlap' in memberships:
        pathways_df = calculate_overlap_membership(pathways_df)
    if 'expansion' in memberships and string_scores_path:
        pathways_df['Expansion_Membership'] = 1
        new_genes_df = calculate_expansion_membership(pathways_df, string_scores_path)
        pathways_df = pd.concat([pathways_df, new_genes_df], ignore_index=True)

    pathways_df.to_csv(output_file, sep="\t", index=False)
    print(f"Filtered and expanded DataFrame saved to {output_file}")

    # Plot histograms if requested
    if plot:
        plot_membership_histograms(pathways_df, os.path.dirname(output_file))
        plot_membership_scatterplots(pathways_df, os.path.dirname(output_file))

    return pathways_df


# Command-line or direct execution
if __name__ == "__main__":
    import sys
    if len(sys.argv) >= 4:
        dataframe_path, output_file = sys.argv[1:3]
        string_scores_path = sys.argv[3] if len(sys.argv) == 5 else None
        memberships = sys.argv[4:-1] if len(sys.argv) > 5 else ['crisp', 'overlap', 'expansion']
        plot = sys.argv[-1].lower() == 'true'
    else:
        dataframe_path = "../../data/pathways/KEGG/KEGG_2022_ensembl.csv" 
        output_file = "../../data/pathways/KEGG/KEGG_2022_membership.tsv" 
        string_scores_path = "../../data/STRING/STRING_GGI.txt" 
        memberships = ['crisp','overlap']
        plot = True

    pathways_df = process_and_expand_pathways(dataframe_path, string_scores_path, output_file, memberships, plot)
