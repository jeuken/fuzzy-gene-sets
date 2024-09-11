#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

def compute_query_membership(adj_p_value):
    """
    Computes the query membership score based on the adjusted p-value using a logistic function.
    
    Args:
        adj_p_value (float): The adjusted p-value to transform.
    
    Returns:
        float: Membership value between 0 and 1, rounded to three decimal places.
    """
    # Calculate the membership score using a logistic transformation
    membership_score = 1 + (-1 / ((1 + np.exp(-50 * adj_p_value)) ** 20))
    return round(membership_score, 3)

def process_query_set(ranked_query_path, output_path):
    """
    Processes the ranked query set to compute Crisp and Fuzzy memberships based on adjusted p-values,
    and saves the results to a specified output path.
    
    Args:
        ranked_query_path (str): Path to the ranked list file containing adjusted p-values.
        output_path (str): Path to save the processed results.
    """
    try:
        # Load the ranked query data with adjusted p-values
        ranked_df = pd.read_csv(ranked_query_path)

        # Compute Crisp Membership: 1 if adj_p_value <= 0.05, otherwise 0
        ranked_df['Crisp_Membership'] = (ranked_df['adj_p_value'] <= 0.05).astype(int)
        
        # Compute Fuzzy Membership values using the membership function
        ranked_df['Fuzzy_Membership'] = ranked_df['adj_p_value'].apply(compute_query_membership)
        
        # Save the results to a file
        ranked_df.to_csv(output_path, sep='\t', index=True)
        print(f"Results saved successfully to {output_path}.")

        # Optionally display the first few rows
        print(ranked_df.head())

    except Exception as e:
        print(f"An error occurred: {e}")

# Application paths
ranked_query_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/gemma/Alzheimers_GSE95587_ranked.csv"
output_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/gemma/Alzheimers_GSE95587_query.csv"

# Process the query set
process_query_set(ranked_query_path, output_path)