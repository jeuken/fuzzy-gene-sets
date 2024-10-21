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

def ora_hypergeometric_test(query_genes_filtered, pathway_genes_filtered, dataset_size):
    """Perform a hypergeometric test to calculate p-values for overlaps."""
    overlap_size = len(set(query_genes_filtered).intersection(set(pathway_genes_filtered)))  # Calculate overlap size
    M = dataset_size  # Total number of genes in the dataset
    n = len(query_genes_filtered)  # Total number of query genes with Membership 1
    N = len(pathway_genes_filtered)  # Total number of genes in the pathway
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
def synthetic_membership(genome_genes, target_genes):
    """Generate a synthetic query set with all target genes and some random genes."""
    genome_genes = list(genome_genes)
    target_genes = list(target_genes)
    
    if len(genome_genes) < 20:  # Ensure enough genes for sampling
        raise ValueError("Not enough genes to sample random genes. Provide a larger genome list.")
    
    # Add all target genes and some random genes from the genome (excluding target genes)
    random_genes = np.random.choice([gene for gene in genome_genes if gene not in target_genes], 
                                    size=int(0.1 * len(target_genes) / (1 - 0.1)), 
                                    replace=False)
    query_genes = set(target_genes).union(set(random_genes))
    
    # Assign membership 1 to genes in query_genes, 0 otherwise
    query_membership = [1 if gene in query_genes else 0 for gene in genome_genes]
    
    return query_membership

# Save synthetic query sets for ORA
def save_synthetic_memberships(n_it, genome_path, pathway_file, pathway_name, synthetic_membership_path):
    """Generate and save synthetic query sets based on specific pathway."""
    # Load genome genes and pathway genes
    genome_genes = pd.read_csv(genome_path, sep='\t')['ENSEMBL'].tolist()
    target_genes = pd.read_csv(pathway_file, sep='\t')
    target_genes = target_genes[target_genes['Pathway_Name'] == pathway_name]['Ensembl_ID'].tolist()
    
    # Generate random memberships
    random_query_memberships_list = []
    for _ in tqdm(range(n_it), desc="Generating synthetic memberships"):
        # Generate random membership for each iteration
        membership = synthetic_membership(genome_genes, target_genes)
        random_query_memberships_list.append(membership)
    
    # Create DataFrame for synthetic memberships
    random_query_memberships = pd.DataFrame(
        data=np.array(random_query_memberships_list).T,  # Transpose to have memberships as columns
        columns=[f"Iter_{i}_Membership" for i in range(n_it)]
    )
    
    # Add Ensembl_ID column
    random_query_memberships['Ensembl_ID'] = genome_genes
    
    # Save synthetic memberships to a file
    random_query_memberships.to_csv(synthetic_membership_path, sep='\t', index=False)
    
    return random_query_memberships

# Run ORA on synthetic query sets to generate null distribution
def standard_ora_null(random_query_memberships_path, pathway_file, n_it, null_distribution_path,pathway_ids):
    """Run ORA on synthetic query sets and collect null distributions."""
    pathways = pd.read_csv(pathway_file, sep='\t')['Pathway_Name'].unique()
    results = pd.DataFrame(index=pathways, columns=[f"Iter_{i}_p_value" for i in range(n_it)])

    def run_iteration(i):
        """Perform standard ORA for one iteration."""
        column_name = f"Iter_{i}_Membership"
        iteration_result = standard_ora(
            query_file=random_query_memberships_path,
            pathway_file=pathway_file,
            query_membership_type=column_name,
            pathway_membership_type='Crisp_Membership',
            output_path=None,
            dataset_name='dataset',
            pathway_ids = pathway_ids
        )
        return iteration_result[['Pathway_Name', 'p_value']]

    # Run iterations in parallel and fill results
    iteration_results = Parallel(n_jobs=-1)(
        delayed(run_iteration)(i) for i in tqdm(range(n_it), desc="Running ORA iterations")
    )

    for i, iteration_result in enumerate(iteration_results):
        column_name = f"Iter_{i}_p_value"
        for pathway_name in iteration_result['Pathway_Name']:
            results.at[pathway_name, column_name] = iteration_result.loc[
                iteration_result['Pathway_Name'] == pathway_name, 'p_value'
            ].values[0]

    # Save null distributions to file
    results.to_csv(null_distribution_path, sep='\t')
    return results

# Parameters and file paths
genome_path = '../../data/IDs/ENSEMBL.txt'
synthetic_membership_path = '../../data/synthetic_query/synthetic_memberships.tsv'
null_distribution_path = '../../data/synthetic_query/null_distributions.csv'
n_it = 10
pathway_file = "../../data/pathways/KEGG/KEGG_2022_membership.tsv"
pathway_name = 'hsa05170'  # Example: Human immunodeficiency virus 1 infection pathway
pathway_ids=  [
    'hsa05166',  # Human T-cell leukemia virus 1 infection
    'hsa05170',  # Human immunodeficiency virus 1 infection
    'hsa05161',  # Hepatitis B
    'hsa05160',  # Hepatitis C
    'hsa05171',  # Coronavirus disease - COVID-19
    'hsa05164',  # Influenza A
    'hsa05162',  # Measles
    'hsa05168',  # Herpes simplex virus 1 infection
    'hsa05163',  # Human cytomegalovirus infection
    'hsa05167',  # Kaposi sarcoma-associated herpesvirus infection
    'hsa05169',  # Epstein-Barr virus infection
    'hsa05165',  # Human papillomavirus infection
    'hsa05110',  # Vibrio cholerae infection
    'hsa05120',  # Epithelial cell signaling in Helicobacter pylori infection
    'hsa05130',  # Pathogenic Escherichia coli infection
    'hsa05132',  # Salmonella infection
    'hsa05131',  # Shigellosis
    'hsa05135',  # Yersinia infection
    'hsa05133',  # Pertussis
    'hsa05134',  # Legionellosis
    'hsa05150',  # Staphylococcus aureus infection
    'hsa05152',  # Tuberculosis
    'hsa05146',  # Amoebiasis
    'hsa05144',  # Malaria
    'hsa05145',  # Toxoplasmosis
    'hsa05140',  # Leishmaniasis
    'hsa05142',  # Chagas disease
    'hsa05143'   # African trypanosomiasis
]
# Step 1: Generate and save synthetic query sets
#save_synthetic_memberships(n_it, genome_path, pathway_file, pathway_name, synthetic_membership_path)

# Step 2: Run ORA null distribution analysis
standard_ora_null(synthetic_membership_path, pathway_file, n_it, null_distribution_path,pathway_ids = pathway_ids)







