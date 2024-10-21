import numpy as np
import pandas as pd
import os
import pandas as pd  # Importing pandas for data manipulation
import scipy.stats as stats  # Importing stats for statistical functions
import os  # Importing os for file and directory operations
import argparse  # Importing argparse for command-line argument parsing
from tqdm import tqdm  # Importing tqdm for progress bars
from joblib import Parallel, delayed,cpu_count

def ora_get_query(query_file, query_membership_type):
    query_df = pd.read_csv(query_file, sep='\t', usecols=['Ensembl_ID', query_membership_type])
    query_df.dropna(subset=[query_membership_type], inplace=True)
    query_df.rename(columns={query_membership_type: 'Membership'}, inplace=True)
    return query_df

def ora_get_pathways(pathway_file, pathway_membership_type):
    df = pd.read_csv(
        pathway_file, 
        sep='\t', 
        usecols=['Pathway_Name', 'Description', 'Ensembl_ID', pathway_membership_type],
        dtype={'Ensembl_ID': str, pathway_membership_type: float}
    ).dropna()
    df.rename(columns={pathway_membership_type: 'Membership'}, inplace=True)
    pathway_df = df.groupby('Pathway_Name').agg({
        'Description': 'first',
        'Ensembl_ID': list,
        'Membership': list
    }).reset_index()
    return pathway_df

def ora_fuzzy_intersection(query_memberships, pathway_memberships):
    intersection = np.multiply(query_memberships, pathway_memberships)  # Element-wise product
    intersection_size = np.sum(intersection)  # Sum of the product values
    return intersection_size


def ora_permutation(query_memberships, pathway_memberships):
    random_intersection_size = ora_fuzzy_intersection(np.random.permutation(query_memberships), pathway_memberships)
    return random_intersection_size

def ora_null_distribution(query_memberships, pathway_memberships, n_permutations=1000, n_jobs=cpu_count()-1):
    null_dist = Parallel(n_jobs=n_jobs)(delayed(ora_permutation)(query_memberships, pathway_memberships) for _ in range(n_permutations))
    return null_dist

def ora_p_value(observed_intersection, null_distribution):
    p_value = np.mean(null_distribution >= observed_intersection)
    return p_value

def ora_compute_stats(pathway, query_df, n_permutations, n_jobs, plots=False):
    pathway_df = pd.DataFrame({
        'Ensembl_ID': pathway['Ensembl_ID'],
        'Membership': pathway['Membership']
    })

    # Perform a left join to keep all query genes and include only matching pathway genes
    merged_df = pd.merge(query_df, pathway_df, on='Ensembl_ID', how='left', suffixes=('_query', '_pathway')).fillna(0)

    # Set pathway memberships to zero for genes not found in the pathway
    merged_df['Membership_pathway'] = merged_df['Membership_pathway'].fillna(0)

    # Get membership arrays for both query and pathway
    query_memberships_full = merged_df['Membership_query'].values
    pathway_memberships_full = merged_df['Membership_pathway'].values

    # Calculate observed intersection
    observed_intersection = ora_fuzzy_intersection(query_memberships_full, pathway_memberships_full)
    
    if plots:
        null_distribution = ora_null_distribution(query_memberships_full, pathway_memberships_full, n_permutations, n_jobs)
        p_value = ora_p_value(observed_intersection, null_distribution)
        return observed_intersection, null_distribution, p_value
    else:
        null_distribution = ora_null_distribution(query_memberships_full, pathway_memberships_full, n_permutations // 10, n_jobs)
        p_value = ora_p_value(observed_intersection, null_distribution)
        return observed_intersection, p_value


def fuzzy_ora(query_file, pathway_file, query_membership_type='Crisp', pathway_membership_type='Crisp', n_permutations=1000, n_jobs=cpu_count()-1, output_path=None, plots=False, dataset_name='dataset', pathway_ids=None):
    query_df = ora_get_query(query_file, query_membership_type)
    pathway_df = ora_get_pathways(pathway_file, pathway_membership_type)
    
    if pathway_ids:
        # Filter pathway_df to include only pathways from the pathway_ids list
        pathway_df = pathway_df[pathway_df['Pathway_Name'].isin(pathway_ids)]
    
    results = []
    num_pathways = len(pathway_df)
    print(f"Total pathways: {num_pathways}")

    for idx, pathway in tqdm(pathway_df.iterrows(), total=num_pathways, desc="Processing Pathways"):
        observed_intersection, p_value = ora_compute_stats(
            pathway, query_df, n_permutations, n_jobs, plots=False)

        results.append({
            'Pathway_Name': pathway['Pathway_Name'],
            'Description': pathway['Description'],
            'Observed_Intersection': observed_intersection,
            'p_value': p_value
        })

    results_df = pd.DataFrame(results)
    results_df.sort_values('p_value', inplace=True)
    results_df.reset_index(drop=True, inplace=True)
    results_df['Rank'] = results_df['p_value'].rank(method='min').astype(int)
    
    if output_path:
        results_folder = os.path.join(output_path, dataset_name)
        os.makedirs(results_folder, exist_ok=True)
        results_file_name = f"ora_results_{dataset_name}_{query_membership_type}_{pathway_membership_type}.csv"
        results_df.to_csv(os.path.join(results_folder, results_file_name), index=False)
    
    return results_df


# Function to generate random memberships
def random_membership(genome_genes):
    dataset_genes = list(genome_genes)
    if len(dataset_genes) < 20:  # Ensure at least 20 genes to sample 5%
        raise ValueError("Not enough genes to sample 5%. Please provide a larger gene list.")
    query_genes = set(np.random.choice(dataset_genes, size=int(0.1 * len(dataset_genes)), replace=False))
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


def fuzzy_ora_null(random_query_memberships_path, pathway_file, n_it, null_distribution_path):
    """Run standard ORA null distribution analysis using random memberships in parallel."""
    pathways = pd.read_csv(pathway_file, sep='\t')['Pathway_Name'].unique()
    results = pd.DataFrame(index=pathways, columns=[f"{i}_p_value" for i in range(n_it)])

    def run_iteration(i):
        """Perform standard ORA for a single iteration and return results."""
        column_name = f"{i}_Membership"
        # Perform standard ORA for the current iteration
        iteration_result = fuzzy_ora(
            query_file=random_query_memberships_path,
            pathway_file=pathway_file,
            query_membership_type=column_name,
            pathway_membership_type='Local_Reaching_Membership',
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
null_distribution_path = '../../data/bias_under_null/null_distributions_crisp_local_reaching.csv'
n_it = 100
pathway_file = "../../data/pathways/KEGG/topology/kegg_topology_memberships.tsv"

# Step 1: Generate and save random memberships
#save_random_memberships(n_it, genome_path, random_membership_path)

# Step 2: Run the standard ORA null distribution analysis with saved memberships
fuzzy_ora_null(random_membership_path, pathway_file, n_it, null_distribution_path)
