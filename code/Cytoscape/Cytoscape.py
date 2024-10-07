#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 11:40:34 2024

@author: piamozdzanowski
"""

import pandas as pd

def csv_to_tsv(input_csv, output_tsv):
    # Read the CSV file
    df = pd.read_csv(input_csv)
    
    # Select the relevant columns
    df = df[["Pathway_Name", "Description", "p_value"]]
    
    # Calculate the rank of the pathway
    df['Rank'] = df['p_value'].rank(method='min',pct = True)
    
    # Calculate p_val as Rank divided by the total number of pathways
    total_pathways = len(df)
    df['p_val'] = df['Rank'] / total_pathways
    
    # Drop the original p_value column
    df.drop(columns='p_value', inplace=True)
    
    # Write the DataFrame to a TSV file (tab-separated)
    df.to_csv(output_tsv, sep='\t', index=False)

# Example usage
csv_to_tsv('/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/ORA_output/T2D_E-HCAD-31/standard_ora_results_T2D_E-HCAD-31_Crisp_Membership_Crisp_Membership.csv', 
            '/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/ORA_output/T2D_E-HCAD-31/standard_ora_results_T2D_E-HCAD-31_Crisp_Membership_Crisp_Membership.tsv')

# Example usage
csv_to_tsv('/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/ORA_output/T2D_E-HCAD-31/ora_results_T2D_E-HCAD-31_Crisp_Membership_Overlap_Membership_linear.csv', 
            '/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/ORA_output/T2D_E-HCAD-31/ora_results_T2D_E-HCAD-31_Crisp_Membership_Overlap_Membership_linear.tsv')
