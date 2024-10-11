import pandas as pd

def create_expression_dataframe(dataset_file, id_mapping_file, output_file):
   
    try:
        id_dict = pd.read_csv(id_mapping_file, sep='\t', header=0, 
                                names=['Ensembl_ID', 'Gene_Name'], dtype=str)
        id_dict.dropna(inplace=True)
        id_dict = id_dict.drop_duplicates(subset='Gene_Name').set_index('Gene_Name')
        print(f"Loaded gene dictionary with {len(id_dict)} entries.")
        
        
        expression_df = pd.read_csv(dataset_file, sep='\t', comment='#')
        initial_row_count = len(expression_df)
        expression_df = expression_df.dropna(subset=['GeneSymbol'])
        dropped_row_count = initial_row_count - len(expression_df)
        print(f"Dropped {dropped_row_count} out of {initial_row_count} rows due to missing GeneSymbols.")
        
        expr_cols = [col for col in expression_df.columns if 'Name=' in col]
        if not expr_cols:
            raise ValueError("No expression columns found in the dataset.")
        
        
        control_keywords = ['healthy', 'con', 'ref', 'wt', 'norm','c_','ctrl']
        

        control_cols = [col for col in expr_cols if any(keyword in col.lower() for keyword in control_keywords)]
        disease_cols = [col for col in expr_cols if col not in control_cols]
        
        if not control_cols or not disease_cols:
            raise ValueError("Not enough columns for control or disease samples.")
        
        expression_df['Ensembl_ID'] = expression_df['GeneSymbol'].map(id_dict['Ensembl_ID'])

        # Drop rows where Ensembl_ID is missing
        initial_row_count = len(expression_df)
        expression_df.dropna(subset=['Ensembl_ID'], inplace=True)
        dropped_row_count = initial_row_count - len(expression_df)
        print(f"Dropped {dropped_row_count} more rows due to missing Ensembl IDs.")
        
        # Group by Ensembl_ID and average the expression values
        expression_df = expression_df.groupby('Ensembl_ID').agg({col: 'mean' for col in control_cols + disease_cols})
        
        # Reset index to make Ensembl_ID a column again
        expression_df.reset_index(inplace=True)
        
        # Rename columns
        control_col_names = [f'control{i+1}' for i in range(len(control_cols))]
        disease_col_names = [f'disease{i+1}' for i in range(len(disease_cols))]
        expression_df.columns = ['Ensembl_ID'] + control_col_names + disease_col_names
        
        # Save the DataFrame to a file
        expression_df.to_csv(output_file, index=False)
        
        print(f"DataFrame saved successfully to {output_file}.")
        
        # Display the first few rows of the DataFrame
        print(expression_df.head())

    except Exception as e:
        print(f"An error occurred: {e}")

# Application
dataset_file = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/gemma/Alzheimers_GSE53697.txt"
id_mapping_file = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/IDs/symbol2ensembl.txt"
output_file = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/gemma/Alzheimers_GSE53697_formatted.csv"

create_expression_dataframe(dataset_file, id_mapping_file, output_file)

