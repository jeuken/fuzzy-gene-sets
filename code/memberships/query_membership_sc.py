import pandas as pd
import numpy as np
import math

def compute_query_membership(q):
    # Fuzzy membership function based on q-value
    fuzzy_membership = 1 - 1 / (1 + np.exp(-300 * (q - 0.05)))
    return fuzzy_membership

def compute_logfc_membership(logfc, logfc_threshold):
    # Fuzzy membership function based on absolute logFC
    fuzzy_membership_logfc = 1 / (1 + np.exp(-300 * (abs(logfc) - logfc_threshold)))
    return fuzzy_membership_logfc

def process_query_set(query_path, output_path, logfc_threshold=1.0):
    # Load the ranked query data with adjusted p-values
    query_df = pd.read_csv(query_path)

    # Compute Crisp Membership with both adjusted p-value and logFC threshold
    query_df['Crisp_Membership'] = ((query_df['q'] <= 0.05) & 
                                      (abs(query_df['logfc']) >= logfc_threshold)).astype(int)

    # Compute Fuzzy Membership values using the membership functions
    query_df['Fuzzy_Membership_Q'] = query_df['q'].apply(compute_query_membership)
    query_df['Fuzzy_Membership_LogFC'] = query_df['logfc'].apply(compute_logfc_membership, logfc_threshold=logfc_threshold)

    # Take the minimum of the two fuzzy memberships
    query_df['Fuzzy_Membership'] = query_df[['Fuzzy_Membership_Q', 'Fuzzy_Membership_LogFC']].min(axis=1)
    
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

    return query_df

# Application paths
query_path = "../../data/single_cell/HIV_E-GEOD-111727.csv"
output_path = "../../data/single_cell/HIV_E-GEOD-111727_membership.csv"

# Process the query set with a specified logFC threshold
logfc_threshold = math.log(1.2)  # Change this value as needed
query_df = process_query_set(query_path, output_path, logfc_threshold)
