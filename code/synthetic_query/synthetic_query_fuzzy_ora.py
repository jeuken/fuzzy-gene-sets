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
def fuzzy_ora_null(random_query_memberships_path, pathway_file, n_it, null_distribution_path,pathway_ids = None):
    """Run ORA on synthetic query sets and collect null distributions."""
    pathways = pd.read_csv(pathway_file, sep='\t')['Pathway_Name'].unique()
    results = pd.DataFrame(index=pathways, columns=[f"Iter_{i}_p_value" for i in range(n_it)])

    def run_iteration(i):
        """Perform standard ORA for one iteration."""
        column_name = f"Iter_{i}_Membership"
        iteration_result = fuzzy_ora(
            query_file=random_query_memberships_path,
            pathway_file=pathway_file,
            query_membership_type=column_name,
            pathway_membership_type='Overlap_Membership',
            output_path=None,
            dataset_name='dataset',
            n_permutations=10000,
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
null_distribution_path = '../../data/synthetic_query/null_distributions_overlap.csv'
n_it = 1
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
fuzzy_ora_null(synthetic_membership_path, pathway_file, n_it, null_distribution_path)







