import pandas as pd
import re

def create_expression_dataframe(dataset_path, mapping_path, output_path):
    """
    Creates a DataFrame with Ensembl IDs, healthy sample columns, and cancer sample columns,
    and saves it to a specified output file.
    
    Args:
        dataset_path (str): Path to the dataset file.
        mapping_path (str): Path to the gene mapping file.
        output_path (str): Path where the DataFrame will be saved.
    """
    try:
        # Load the gene mapping file with appropriate header handling
        gene_dict = pd.read_csv(mapping_path, sep='\t', header=0, 
                                names=['Gene_ID', 'Gene_Name'], dtype=str)
        
        # Drop rows with NaN values and duplicates
        gene_dict.dropna(inplace=True)
        gene_dict.drop_duplicates(subset='Gene_Name', inplace=True)
        gene_dict.set_index('Gene_Name', inplace=True)
        
        print(f"Loaded gene dictionary with {len(gene_dict)} entries.")
        
        # Load the dataset, skipping lines starting with '#'
        df = pd.read_csv(dataset_path, sep='\t', comment='#')
        
        # Drop rows where GeneSymbol is missing
        initial_row_count = len(df)
        df = df.dropna(subset=['GeneSymbol'])
        dropped_row_count = initial_row_count - len(df)
        print(f"Rows dropped due to missing GeneSymbols: {dropped_row_count}")
        
    
        # Identify columns containing expression data
        expr_cols = [col for col in df.columns if 'Name=' in col]
        if not expr_cols:
            raise ValueError("No expression columns found in the dataset.")
        
        # Extract healthy and cancer samples
        healthy_cols = [col for col in expr_cols if 'Normal' in col]
        cancer_cols = [col for col in expr_cols if 'RenalClearCellCarcinoma' in col]
        
        if not healthy_cols or not cancer_cols:
            raise ValueError("Not enough columns for healthy or cancer samples.")
        
        # Map gene symbols to Ensembl IDs and create the new DataFrame
        df['Ensembl_ID'] = df['GeneSymbol'].map(gene_dict['Gene_ID'])
        
        # Drop rows where Ensembl_ID is missing
        initial_row_count = len(df)
        df.dropna(subset=['Ensembl_ID'], inplace=True)
        dropped_row_count = initial_row_count - len(df)
        print(f"Rows dropped due to missing Ensembl IDs: {dropped_row_count}")
        
        # Group by Ensembl_ID and average the expression values
        df = df.groupby('Ensembl_ID').agg({col: 'mean' for col in healthy_cols + cancer_cols})
        
        # Reset index to make Ensembl_ID a column again
        df.reset_index(inplace=True)
        
        # Rename columns to match the required format
        healthy_col_names = [f'healthy{i+1}' for i in range(len(healthy_cols))]
        cancer_col_names = [f'cancer{i+1}' for i in range(len(cancer_cols))]
        df.columns = ['Ensembl_ID'] + healthy_col_names + cancer_col_names
        
        # Save the DataFrame to a file
        df.to_csv(output_path, sep='\t', index=False)
        
        print(f"DataFrame saved successfully to {output_path}.")
        
        # Display the first few rows of the DataFrame
        print(df.head())
        
    except Exception as e:
        print(f"An error occurred: {e}")

# Example usage
dataset_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/gemma/practice_renal_carcinoma_samples.txt"
mapping_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/IDs/symbol2ensembl.txt"
output_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/gemma/expression_dataframe.txt"

create_expression_dataframe(dataset_path, mapping_path, output_path)
