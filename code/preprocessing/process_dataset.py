import pandas as pd

def create_expression_dataframe(dataset_path, mapping_path, output_path):
    """
    Creates a DataFrame with Ensembl IDs, control sample columns, and disease sample columns,
    and saves it to a specified output file.

    Args:
        dataset_path (str): Path to the dataset file.
        mapping_path (str): Path to the gene mapping file.
        output_path (str): Path where the DataFrame will be saved.
    """
    try:
        # Load the gene mapping file
        gene_dict = pd.read_csv(mapping_path, sep='\t', header=0, 
                                names=['Ensembl_ID', 'Gene_Name'], dtype=str)
        
        # Drop rows with NaN values and ensure unique Gene_Name
        gene_dict.dropna(inplace=True)
        gene_dict = gene_dict.drop_duplicates(subset='Gene_Name').set_index('Gene_Name')

        print(f"Loaded gene dictionary with {len(gene_dict)} entries.")
        
        # Load the dataset, skipping lines starting with '#'
        df = pd.read_csv(dataset_path, sep='\t', comment='#')
        
        # Drop rows where GeneSymbol is missing
        initial_row_count = len(df)
        df = df.dropna(subset=['GeneSymbol'])
        dropped_row_count = initial_row_count - len(df)
        print(f"Dropped {dropped_row_count} rows due to missing GeneSymbols.")
        
        # Identify columns containing expression data
        expr_cols = [col for col in df.columns if 'Name=' in col]
        if not expr_cols:
            raise ValueError("No expression columns found in the dataset.")
        
        # Define keywords for control samples
        control_keywords = ['healthy', 'con', 'ref', 'wt', 'norm']
        
        # Assign columns to control or disease based on keywords
        control_cols = [col for col in expr_cols if any(keyword in col.lower() for keyword in control_keywords)]
        disease_cols = [col for col in expr_cols if col not in control_cols]
        
        if not control_cols or not disease_cols:
            raise ValueError("Not enough columns for control or disease samples.")
        
        # Map gene symbols to Ensembl IDs and create the new DataFrame
        df['Ensembl_ID'] = df['GeneSymbol'].map(gene_dict['Ensembl_ID'])

        # Drop rows where Ensembl_ID is missing
        initial_row_count = len(df)
        df.dropna(subset=['Ensembl_ID'], inplace=True)
        dropped_row_count = initial_row_count - len(df)
        print(f"Dropped {dropped_row_count} rows due to missing Ensembl IDs.")
        
        # Group by Ensembl_ID and average the expression values
        df = df.groupby('Ensembl_ID').agg({col: 'mean' for col in control_cols + disease_cols})
        
        # Reset index to make Ensembl_ID a column again
        df.reset_index(inplace=True)
        
        # Rename columns
        control_col_names = [f'control{i+1}' for i in range(len(control_cols))]
        disease_col_names = [f'disease{i+1}' for i in range(len(disease_cols))]
        df.columns = ['Ensembl_ID'] + control_col_names + disease_col_names
        
        # Save the DataFrame to a file
        df.to_csv(output_path, index=False)
        
        print(f"DataFrame saved successfully to {output_path}.")
        
        # Display the first few rows of the DataFrame
        print(df.head())

    except Exception as e:
        print(f"An error occurred: {e}")

# Application
dataset_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/gemma/Alzheimers_GSE95587.txt"
mapping_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/IDs/symbol2ensembl.txt"
output_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/gemma/Alzheimers_GSE95587_formatted.csv"

create_expression_dataframe(dataset_path, mapping_path, output_path)

