import pandas as pd
import numpy as np

def compute_query_membership(q):
    # Fuzzy membership function
    fuzzy_membership = 1 - 1 / (1 + np.exp(-300 * (q - 0.05)))  # Corrected multiplication with *
    return fuzzy_membership

def process_query_set(query_path, output_path):
    # Load the ranked query data with adjusted p-values
    query_df = pd.read_csv(query_path)
    
    # Compute Crisp Membership
    query_df['Crisp_Membership'] = (query_df['adj_p_value'] <= 0.05).astype(int)

    # Compute Fuzzy Membership values using the membership function
    query_df['Fuzzy_Membership'] = query_df['adj_p_value'].apply(compute_query_membership)
    
    # Save the results to a file
    query_df.to_csv(output_path, sep='\t', index=False)
    print(f"Results saved successfully to {output_path}.")

    # Calculate the number and percentage of genes with Crisp Membership of 1
    num_crisp_membership_1 = query_df['Crisp_Membership'].sum()
    total_genes = len(query_df)
    percentage_crisp_membership_1 = (num_crisp_membership_1 / total_genes) * 100 if total_genes > 0 else 0

    # Print the results
    print(f"Number of genes with Crisp Membership of 1: {num_crisp_membership_1}")
    print(f"Percentage of genes with Crisp Membership of 1: {percentage_crisp_membership_1:.2f}%")

    # Optionally display the first few rows
    print(query_df.head())
    return query_df

# Application paths
query_path = "../../data/cancer/breast_cancer/breast_cancer_GSE9574_t_stats.csv"
output_path = "../../data/cancer/breast_cancer/breast_cancer_GSE9574_membership.csv"

# Process the query set
query_df = process_query_set(query_path, output_path)
