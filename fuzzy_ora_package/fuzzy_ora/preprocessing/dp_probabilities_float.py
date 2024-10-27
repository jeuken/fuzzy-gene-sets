import os
import numpy as np
import pandas as pd
from joblib import Parallel, delayed, cpu_count
from tqdm import tqdm
from scipy.special import comb

def genesetdp(k, step_size=1/337):
    max_s = sum(m * c for m, c in k)  
    Q = sum(c for _, c in k)  
    scaled_max_s = int(max_s / step_size)  

    N = np.zeros((scaled_max_s + 1, Q + 1), dtype=np.float64)  
    N[0, 0] = 1  

    min_membership_scaled = int(min(m for m, _ in k) / step_size)

    for membership_value, count in k:
        scaled_membership = int(membership_value / step_size)  
        
        for s in range(scaled_max_s, scaled_membership - 1, -1):  
            min_genes = max(0, (s + scaled_membership - 1) // scaled_membership)  
            max_genes = min(Q, s // min_membership_scaled)  
            
            for c in range(min_genes, max_genes + 1):  
                for b in range(1, min(count, c) + 1):  
                    if s >= b * scaled_membership:
                        N[s, c] += comb(count, b,exact = True) * N[s - b * scaled_membership, c - b]

    row_index = [i * step_size for i in range(scaled_max_s + 1)]
    return pd.DataFrame(N, index=row_index, columns=range(Q + 1))


# Calculate probabilities of observing a score >= s given c
def calculate_probabilities(N_df):
    prob_df = pd.DataFrame(index=N_df.index, columns=N_df.columns)
    
    for c in N_df.columns:
        cumulative_counts = N_df[c].cumsum().values
        total_ways = cumulative_counts[-1]
        prob_array = np.zeros_like(cumulative_counts)
        
        for s in range(len(prob_array)):
            if total_ways > 0:
                prob_array[s] = (total_ways - cumulative_counts[s-1]) / total_ways if s > 0 else total_ways / total_ways
            else:
                prob_array[s] = 0

        prob_df[c] = prob_array

    return prob_df

# Function to create membership tuples for each pathway
def create_membership_tuples(df):
    pathway_membership_dict = {}
    for pathway, group in df.groupby('Pathway_Name'):
        membership_counts = group.groupby('Overlap_Membership')['Ensembl_ID'].count().reset_index()
        membership_tuples = list(membership_counts.itertuples(index=False, name=None))
        membership_tuples.sort(key=lambda x: x[0])
        pathway_membership_dict[pathway] = membership_tuples
    return pathway_membership_dict

# Save results for each pathway (parallelized)
def process_pathway(pathway, k, output_folder):
    N_df = genesetdp(k, step_size=1/337)  
    probabilities_df = calculate_probabilities(N_df)
    
    output_file_path = os.path.join(output_folder, f'{pathway}_dp_probabilities.tsv')
    
    probabilities_df.to_csv(output_file_path, sep='\t')
    
    return f"Saved probabilities for pathway: {pathway} to {output_file_path}"

# Main function
def main(selected_pathways=None):
    # Load the TSV file into a DataFrame
    df = pd.read_csv('../../../data/pathways/KEGG/KEGG_2022_pathway_memberships.tsv', sep='\t')

    # Generate the dictionary of pathway membership tuples
    pathway_memberships = create_membership_tuples(df)

    # Filter pathways if selected_pathways is provided
    if selected_pathways is not None:
        pathway_memberships = {pathway: k for pathway, k in pathway_memberships.items() if pathway in selected_pathways}

    # Create an output folder to save the results if it doesn't exist
    output_folder = '../../../data/pathways/KEGG/DP/Overlap_2/'
    os.makedirs(output_folder, exist_ok=True)

    # Use joblib for parallel processing and tqdm for the progress bar
    Parallel(n_jobs = cpu_count()-1)(
        delayed(process_pathway)(pathway, k, output_folder)
        for pathway, k in tqdm(pathway_memberships.items(), desc="Processing pathways")
    )
    
if __name__ == "__main__":
    # Example pathways to process
    specified_pathways = None  # Customize this list
    main(selected_pathways=specified_pathways)

