import numpy as np
import pandas as pd
import os
import pandas as pd  # Importing pandas for data manipulation
import scipy.stats as stats  # Importing stats for statistical functions
import os  # Importing os for file and directory operations
import argparse  # Importing argparse for command-line argument parsing
from tqdm import tqdm  # Importing tqdm for progress bars
from joblib import Parallel, delayed

def ora_get_query(query_file, query_membership_type):
    """Load the query file and filter it based on the specified membership type."""
    if not os.path.isfile(query_file):
        raise FileNotFoundError(f"The query file {query_file} does not exist.")  # Check if the query file exists
    query_df = pd.read_csv(query_file, sep='\t', usecols=['Ensembl_ID', query_membership_type])  # Read the query file
    query_df.dropna(subset=[query_membership_type], inplace=True)  # Drop rows with NaN values in Membership
    query_df.rename(columns={query_membership_type: 'Membership'}, inplace=True)  # Rename the membership column
    return query_df  # Return the filtered DataFrame

def ora_get_pathways(pathway_file, pathway_membership_type):
    """Load the pathway file and group pathways by their names."""
    df = pd.read_csv(
        pathway_file, 
        sep='\t', 
        usecols=['Pathway_Name', 'Description', 'Ensembl_ID', pathway_membership_type],
        dtype={'Ensembl_ID': str, pathway_membership_type: float}
    ).dropna()  # Read and drop NaN values
    df.rename(columns={pathway_membership_type: 'Membership'}, inplace=True)  # Rename the membership column
    # Group by Pathway_Name and aggregate the data
    pathway_df = df.groupby('Pathway_Name').agg({
        'Description': 'first',
        'Ensembl_ID': list,
        'Membership': list
    }).reset_index()  # Reset index after grouping
    return pathway_df  # Return the aggregated DataFrame

def ora_hypergeometric_test(query_genes, pathway_genes, dataset_size):
    """Perform a hypergeometric test to calculate p-values for overlaps."""
    overlap_size = len(set(query_genes).intersection(set(pathway_genes)))  # Calculate overlap size
    M = dataset_size  # Total number of genes in the dataset
    n = len(query_genes)  # Total number of query genes with Membership 1
    N = len(pathway_genes)  # Total number of genes in the pathway
    k = overlap_size  # Number of overlapping genes

    # Calculate p-value using the hypergeometric distribution
    p_value = stats.hypergeom.sf(k - 1, M, N, n)
    
    return overlap_size, p_value  # Return overlap size and p-value

def standard_ora_compute_stats(pathway, query_df, dataset_size):
    """Compute statistics for the standard ORA analysis."""
    query_genes = query_df['Ensembl_ID'].tolist()  # List of all query genes
    query_genes_filtered = query_df[query_df['Membership'] == 1]['Ensembl_ID'].tolist()  # Filter query genes by Membership

    pathway_genes = pathway['Ensembl_ID']  # List of Ensembl IDs in the pathway
    pathway_memberships = pathway['Membership']  # List of Memberships in the pathway

    # Filter pathway genes based on Membership and presence in query_genes
    pathway_genes_filtered = [
        gene for gene, membership in zip(pathway_genes, pathway_memberships) 
        if membership == 1 and gene in query_genes  # Must have Membership 1 and be in query_genes
    ]

    overlap_size, p_value = ora_hypergeometric_test(query_genes_filtered, pathway_genes_filtered, dataset_size)  # Perform the hypergeometric test
    
    return overlap_size, p_value  # Return overlap size and p-value

def standard_ora(query_file, pathway_file, query_membership_type='Crisp_Membership', pathway_membership_type='Crisp_Membership', output_path=None, dataset_name='dataset', pathway_ids=None):
    """Main function to perform standard ORA analysis."""
    query_df = ora_get_query(query_file, query_membership_type)  # Load and filter query data
    pathway_df = ora_get_pathways(pathway_file, pathway_membership_type)  # Load and group pathway data
    
    # Filter pathways if specific pathway IDs are provided
    if pathway_ids:
        pathway_df = pathway_df[pathway_df['Pathway_Name'].isin(pathway_ids)]
    
    dataset_size = len(query_df)  # Total number of genes in the query file

    results = []  # Initialize a list to store results
    num_pathways = len(pathway_df)  # Total number of pathways

    # Iterate through each pathway and compute stats
    for idx, pathway in pathway_df.iterrows():
        pathway_genes = pathway['Ensembl_ID']  # Get pathway genes
        
        # Filter pathway genes that are present in the query
        pathway_genes_filtered = [gene for gene in pathway_genes if gene in query_df['Ensembl_ID'].values]
        
        # Only consider pathways with more than 15 genes
        if len(pathway_genes_filtered) > 15:
            overlap_size, p_value = standard_ora_compute_stats(pathway, query_df, dataset_size)  # Compute statistics
            
            # Append results as a dictionary
            results.append({
                'Pathway_Name': pathway['Pathway_Name'],
                'Description': pathway['Description'],
                'Observed_Intersection': overlap_size,
                'p_value': p_value
            })

    results_df = pd.DataFrame(results)  # Convert results to DataFrame
    results_df.sort_values('p_value', inplace=True)  # Sort by p-value
    results_df.reset_index(drop=True, inplace=True)  # Reset index
    results_df['Rank'] = results_df['p_value'].rank(method='min').astype(int)  # Assign ranks based on p-value

    # Save results to CSV if output path is provided
    if output_path:
        results_folder = os.path.join(output_path, dataset_name)
        os.makedirs(results_folder, exist_ok=True)  # Create output directory if it doesn't exist
        results_file_name = f"standard_ora_results_{dataset_name}_{query_membership_type}_{pathway_membership_type}.csv"
        results_df.to_csv(os.path.join(results_folder, results_file_name), index=False)  # Save results to CSV
    
    return results_df  # Return the results DataFrame


# Function to generate random memberships
def random_membership(genome_genes):
    dataset_genes = list(genome_genes)
    if len(dataset_genes) < 20:  # Ensure at least 20 genes to sample 5%
        raise ValueError("Not enough genes to sample 5%. Please provide a larger gene list.")
    query_genes = set(np.random.choice(dataset_genes, size=int(0.3 * len(dataset_genes)), replace=False))
    random_query_membership = [1 if gene in query_genes else 0 for gene in genome_genes]
    return random_query_membership

def save_random_memberships(n_it, genome_path, random_membership_path):
    # Read genes from file and convert to a list
    genome_genes = pd.read_csv(genome_path, sep='\t')['ENSEMBL'].tolist()

    # Generate random memberships
    random_query_memberships_list = []
    for _ in tqdm(range(n_it), desc="Generating random memberships"):
        # Generate a random membership for each iteration
        membership = random_membership(genome_genes)
        random_query_memberships_list.append(membership)

    # Create a DataFrame with the random memberships
    random_query_memberships = pd.DataFrame(
        data=np.array(random_query_memberships_list).T,  # Transpose the array
        columns=[f"{i}_Membership" for i in range(n_it)]
    )

    # Add Ensembl_ID as a column
    random_query_memberships['Ensembl_ID'] = genome_genes

    # Save the random memberships to a CSV file
    random_query_memberships.to_csv(random_membership_path, sep='\t', index=False)
    
    return random_query_memberships


def standard_ora_null(random_query_memberships_path, pathway_file, n_it, null_distribution_path):
    """Run standard ORA null distribution analysis using random memberships in parallel."""
    pathways = pd.read_csv(pathway_file, sep='\t')['Pathway_Name'].unique()
    results = pd.DataFrame(index=pathways, columns=[f"{i}_p_value" for i in range(n_it)])

    def run_iteration(i):
        """Perform standard ORA for a single iteration and return results."""
        column_name = f"{i}_Membership"
        # Perform standard ORA for the current iteration
        iteration_result = standard_ora(
            query_file=random_query_memberships_path,
            pathway_file=pathway_file,
            query_membership_type=column_name,
            pathway_membership_type='Crisp_Membership',
            output_path=None,
            dataset_name='dataset',
            pathway_ids=None
        )
        
        # Collect results for this iteration
        return iteration_result[['Pathway_Name', 'p_value']]

    # Run iterations in parallel with a progress bar
    iteration_results = Parallel(n_jobs=-1)(delayed(run_iteration)(i) for i in tqdm(range(n_it), desc="Running ORA iterations"))

    # Fill the results DataFrame with p-values
    for i, iteration_result in enumerate(iteration_results):
        column_name = f"{i}_p_value"
        for pathway_name in iteration_result['Pathway_Name']:
            results.at[pathway_name, column_name] = iteration_result.loc[iteration_result['Pathway_Name'] == pathway_name, 'p_value'].values[0]

    # Save the results to a file
    results.to_csv(null_distribution_path, sep='\t')
    return results

# Paths and parameters
genome_path = '../../data/IDs/ENSEMBL.txt'
#genome_path = '/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/IDs/ENTREZID_to_ENSEMBL.txt'

random_membership_path = '../../data/bias_under_null/random_memberships.tsv'
null_distribution_path = '../../data/bias_under_null/null_distributions.csv'
n_it = 10
pathway_file = "../../data/pathways/KEGG/KEGG_2022_membership.tsv"

# Step 1: Generate and save random memberships
save_random_memberships(n_it, genome_path, random_membership_path)

# Step 2: Run the standard ORA null distribution analysis with saved memberships
standard_ora_null(random_membership_path, pathway_file, n_it, null_distribution_path)
