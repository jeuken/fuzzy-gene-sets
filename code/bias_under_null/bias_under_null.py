import numpy as np
import pandas as pd
from tqdm import tqdm
import os

# Set the working directory
os.chdir('/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/code/')
from ora.standard_ora_new import standard_ora

# Function to generate random memberships
def random_membership(genome_genes):
    dataset_genes = list(genome_genes)
    if len(dataset_genes) < 20:  # Ensure at least 20 genes to sample 5%
        raise ValueError("Not enough genes to sample 5%. Please provide a larger gene list.")
    query_genes = set(np.random.choice(dataset_genes, size=int(0.1 * len(dataset_genes)), replace=False))
    random_query_membership = [1 if gene in query_genes else 0 for gene in genome_genes]
    return random_query_membership

# Function to save random memberships
def save_random_memberships(n_it, genome_path, random_membership_path):
    # Read genes from file and convert to a list
    genome_genes = pd.read_csv(genome_path, header=None)[0].tolist()

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
    # Get unique pathway names
    pathways = pd.read_csv(pathway_file, sep='\t')['Pathway_Name'].unique()
    results = pd.DataFrame(index=pathways, columns=[f"{i}_p_value" for i in range(n_it)])

    for i in tqdm(range(n_it), desc="Running ORA iterations"):
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
        
        # Check if 'p_value' column exists in the iteration result
        if 'p_value' in iteration_result.columns:
            # Use 'Pathway_Name' to ensure p-values are assigned correctly
            for pathway_name in iteration_result['Pathway_Name']:
                results.at[pathway_name, f"{i}_p_value"] = iteration_result.loc[iteration_result['Pathway_Name'] == pathway_name, 'p_value'].values[0]
        else:
            print(f"Warning: 'p_value' not found in iteration {i} result.")

    # Save the results to a file
    results.to_csv(null_distribution_path, sep='\t')
    return results


# Paths and parameters
genome_path = '../data/IDs/ENSEMBL.txt'
random_membership_path = '../data/bias_under_null/random_memberships.tsv'
null_distribution_path = '../data/bias_under_null/null_distributions.csv'
n_it = 1000
pathway_file = "../data/pathways/KEGG/KEGG_2022_membership.tsv"

# Step 1: Generate and save random memberships
save_random_memberships(n_it, genome_path, random_membership_path)

# Step 2: Run the standard ORA null distribution analysis with saved memberships
standard_ora_null(random_membership_path, pathway_file, n_it, null_distribution_path)
