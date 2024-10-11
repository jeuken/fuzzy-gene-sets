import pandas as pd
import numpy as np
import math

def compute_query_membership(q):
    fuzzy_membership = 1 + (-1 / ((1 + np.exp(-50 * q)) ** 20))
    return round(fuzzy_membership, 3)

def process_query_set(query_path, output_path):
    # Load the ranked query data with adjusted p-values
    query_df = pd.read_csv(query_path)
    
    # Define your logFC threshold
    logfc_threshold = math.log(1.1) 
    
    # Compute Crisp Membership
    query_df['Crisp_Membership'] = ((query_df['q'] <= 0.05) & (query_df['logfc'].abs() >= logfc_threshold)).astype(int)

    # Compute Fuzzy Membership values using the membership function
    query_df['Fuzzy_Membership'] = query_df['q'].apply(compute_query_membership)
    
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
query_path = "../../data/single_cell/T2D_E-HCAD-31.csv"
output_path = "../../data/single_cell/T2D_E-HCAD-31_membership.csv"

# Process the query set
query_df = process_query_set(query_path, output_path)
