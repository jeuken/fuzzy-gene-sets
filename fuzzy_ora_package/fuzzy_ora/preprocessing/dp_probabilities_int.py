import os
import numpy as np
import pandas as pd
from joblib import Parallel, delayed, cpu_count
from tqdm import tqdm
from scipy.special import comb
from decimal import Decimal, getcontext

def genesetdp(k):
    max_s = sum(m * c for m, c in k)  # Sum of membership * count for maximum possible score
    Q = sum(c for _, c in k)          # Total count of genes

    # Use np.object_ to store large integers
    N = np.zeros((max_s + 1, Q + 1), dtype=np.object_)
    N[0, 0] = 1  

    min_membership = min(m for m, _ in k)  # Minimum membership value

    for membership_value, count in k:
        for s in range(max_s, membership_value - 1, -1):
            min_genes = max(0, (s + membership_value - 1) // membership_value)  # Minimum gene count that can achieve `s`
            max_genes = min(Q, s // min_membership)  # Maximum gene count that can achieve `s`
            
            for c in range(min_genes, max_genes + 1):
                for b in range(1, min(count, c) + 1):
                    if s >= b * membership_value:
                        N[s, c] += comb(count, b, exact=True) * N[s - b * membership_value, c - b]

    row_index = list(range(max_s + 1))  # Integer indices, since we’re not scaling
    return pd.DataFrame(N, index=row_index, columns=range(Q + 1))


# Calculate probabilities of observing a score >= s given c
def calculate_probabilities(N_df):
    prob_df = pd.DataFrame(index=N_df.index, columns=N_df.columns, dtype=object)
    
    for c in N_df.columns:
        cumulative_counts = N_df[c].cumsum().astype(object)  # Keeps large integers in object format
        total_ways = Decimal(cumulative_counts.iloc[-1])  # Convert to Decimal for precise division
        
        prob_array = np.zeros_like(cumulative_counts, dtype=object)
        
        for s in range(len(prob_array)):
            if total_ways > 0:
                prob_array[s] = (total_ways - Decimal(cumulative_counts.iloc[s-1])) / total_ways if s > 0 else Decimal(1)
            else:
                prob_array[s] = Decimal(0)

        prob_df[c] = prob_array
    return prob_df

def create_membership_tuples(df,factor = 337):
    pathway_membership_dict = {}
    
    for pathway, group in df.groupby('Pathway_Name'):
        # Scale memberships by fmax and convert to integers
        group['Scaled_Membership'] = (group['Overlap_Membership'] * factor).round().astype(int)
        
        # Count occurrences for each scaled membership value
        membership_counts = group.groupby('Scaled_Membership')['Ensembl_ID'].count().reset_index()
        
        # Create tuples and sort by scaled membership
        membership_tuples = list(membership_counts.itertuples(index=False, name=None))
        membership_tuples.sort(key=lambda x: x[0])  # Sort by scaled membership value
        pathway_membership_dict[pathway] = membership_tuples
    
    return pathway_membership_dict


# Save results for each pathway (parallelized)
def process_pathway(pathway, k, output_folder):
    N_df = genesetdp(k)  
    probabilities_df = calculate_probabilities(N_df)
    
    output_file_path = os.path.join(output_folder, f'{pathway}_dp_probabilities.tsv')
    probabilities_df.to_csv(output_file_path, sep='\t')
    
    return f"Saved probabilities for pathway: {pathway} to {output_file_path}"

# Main function
def main(pathway_file,selected_pathways=None,factor = 337):
    # Load the TSV file into a DataFrame
    df = pd.read_csv(pathway_file, sep='\t')

    # Generate the dictionary of pathway membership tuples
    pathway_memberships = create_membership_tuples(df,factor)

    # Filter pathways if selected_pathways is provided
    if selected_pathways is not None:
        pathway_memberships = {pathway: k for pathway, k in pathway_memberships.items() if pathway in selected_pathways}

    # Create an output folder to save the results if it doesn't exist
    output_folder = '../../../data/pathways/KEGG/DP/Overlap_int/'
    os.makedirs(output_folder, exist_ok=True)

    # Use joblib for parallel processing and tqdm for the progress bar
    Parallel(n_jobs=cpu_count()-1)(
        delayed(process_pathway)(pathway, k, output_folder)
        for pathway, k in tqdm(pathway_memberships.items(), desc="Processing pathways")
    )
    
if __name__ == "__main__":
    pathway_file = '../../../data/pathways/KEGG/KEGG_2022_pathway_memberships.tsv'
    specified_pathways = None
    factor = 337
    main(pathway_file, selected_pathways=specified_pathways,factor = factor)
