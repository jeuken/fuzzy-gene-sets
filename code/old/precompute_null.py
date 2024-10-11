import numpy as np
import pandas as pd
from joblib import Parallel, delayed, cpu_count
import os
import tqdm
import pickle

# Function to load pathways and memberships from file
def get_pathways(pathway_file, pathway_membership_type):
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

# Function to compute null distributions for a subset of pathways and save them individually
def get_nulls(pathway_file, pathway_membership_type, subset_pathways, n_jobs, n_its, output_dir):
    # Load pathways and memberships
    pathway_df = get_pathways(pathway_file, pathway_membership_type)
    
    # Filter by subset of pathways if provided
    if subset_pathways is not None:
        pathway_df = pathway_df[pathway_df['Pathway_Name'].isin(subset_pathways)]
    
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Loop through each pathway
    for idx, pathway in enumerate(tqdm.tqdm(pathway_df['Pathway_Name'], desc="Processing Pathways")):
        memberships = pathway_df.loc[pathway_df["Pathway_Name"] == pathway, 'Membership'].values[0]
        null_dist = {}

        # For each k from 1 to the length of memberships, compute the null distribution
        for k in range(1, len(memberships) + 1):
            null_dist[k] = compute_null(memberships, k, n_jobs, n_its)

        return null_dist
        # Save the null distribution to a pickle file
        output_path = os.path.join(output_dir, f"{pathway}.pkl")
        with open(output_path, 'wb') as f:
            pickle.dump(null_dist, f)

# Function to compute the null distribution for a specific k
def compute_null(memberships, k, n_jobs, n_its):
    # Use parallel processing to generate the null distribution
    null_dist = Parallel(n_jobs=n_jobs)(
        delayed(permutation)(k, memberships) for _ in range(n_its)
    )
    return null_dist

# Function to perform a single permutation and compute the fuzzy score for k
def permutation(k, memberships):
    # Randomly shuffle the memberships
    permuted_memberships = np.random.permutation(memberships)
    
    # Take the top k elements and sum them (fuzzy intersection)
    score = np.sum(permuted_memberships[:k])
    return score

# Example usage
if __name__ == "__main__":
    pathway_file = "../../data/pathways/KEGG/KEGG_2022_membership.tsv"
    pathway_membership_type = 'Overlap_Membership'
    n_jobs = cpu_count() - 1
    n_its = 100000
    subset_pathways = ['hsa05110','hsa05161']  # Example subset
    output_dir = "../../data/pathways/KEGG/Null_Distributions/Membership_Type"  # Path to save null distributions

    # Compute and save null distributions for the subset of pathways
    null_dist = get_nulls(pathway_file, pathway_membership_type, subset_pathways, n_jobs, n_its, output_dir)



'''import numpy as np
from collections import Counter

def compute_null(memberships, k):
    # Create a counter for the unique membership values
    membership_counter = Counter(memberships)
    
    # Extract unique values and their frequencies
    unique_memberships = list(membership_counter.keys())
    frequencies = list(membership_counter.values())

    # Initialize a dictionary to count occurrences of each sum
    sum_distribution = Counter()

    # Use combinations with repetition to find all sums of k memberships
    def find_sums(current_sum, depth):
        if depth == k:
            sum_distribution[current_sum] += 1
            return
        for value, freq in zip(unique_memberships, frequencies):
            if freq > 0:  # Ensure there are still frequencies available
                find_sums(current_sum + value, depth + 1)

    # Start finding sums
    find_sums(0, 0)

    # Convert Counter to a list of sums and their occurrences
    null_dist = dict(sum_distribution)
    return null_dist

# Example usage
memberships = [0.1, 0.2, 0.3] * 100  # Example membership values
k = 5  # Choosing 5 memberships
null_distribution = compute_null(memberships, k)
print(null_distribution)
'''