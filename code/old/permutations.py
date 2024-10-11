#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 15:47:40 2024

@author: piamozdzanowski
"""

import pandas as pd

# Load the data from the file (assuming the file is named 'expression_data.txt')
file_path = '/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/gemma/practice_renal_carcinoma_samples.txt'  # Replace with the actual file path
df = pd.read_csv(file_path, sep='\t', comment='#')  # Assuming tab-separated values and header lines start with '#'

# Extract 'GeneSymbol' and expression value columns
expression_columns = [col for col in df.columns if 'GSE781_Biomat' in col]
df_expression = df[['GeneSymbol'] + expression_columns]

# Rename expression columns to 'Healthy_X' and 'Cancer_X'
healthy_cols = [col for col in expression_columns if 'NormalHumanKidney' in col]
cancer_cols = [col for col in expression_columns if 'RenalClearCellCarcinoma' in col]

# Renaming columns sequentially
df_expression.rename(columns={col: f'Healthy_{i+1}' for i, col in enumerate(healthy_cols)}, inplace=True)
df_expression.rename(columns={col: f'Cancer_{i+1}' for i, col in enumerate(cancer_cols)}, inplace=True)

# Display the first few rows of the DataFrame
print(df_expression.head())
