import pandas as pd
import os

def format_expression_data(dataset_file, id_mapping_file, output_file):
    # Extract ID types from the filename
    base_filename = os.path.basename(id_mapping_file)
    id_types = base_filename.replace('.txt', '').split('_to_')  # Assuming format is "ID1_to_ID2.txt"
    
    if len(id_types) != 2:
        raise ValueError("ID mapping filename must be in the format 'ID1_to_ID2.txt'")
    
    input_id_type = id_types[0]
    output_id_type = id_types[1]
    
    # Load gene mapping file, drop rows with NA values
    id_dict = pd.read_csv(id_mapping_file, sep='\t', header=0, dtype=str)
    
    # Drop rows where either input_id_type or output_id_type is NA
    id_dict = id_dict.dropna(subset=[input_id_type, output_id_type])
    
    # Keep the first Ensembl ID for duplicate Gene Names
    id_dict = id_dict.drop_duplicates(subset=[input_id_type], keep='first')

    # Rename columns based on extracted ID types
    id_dict.columns = [input_id_type, output_id_type]
    
    # Load dataset and record initial row count
    expression_df = pd.read_csv(dataset_file, sep='\t', comment='#')
    initial_rows = len(expression_df)

    # Drop rows with missing GeneSymbol
    expression_df = expression_df.dropna(subset=['GeneSymbol'])
    gene_symbol_filtered_rows = len(expression_df)

    # Identify expression columns based on keyword patterns
    expr_cols = [col for col in expression_df.columns if 'Name=' in col]
    if not expr_cols:
        raise ValueError("No expression columns found in the dataset.")

    # Define control keywords and separate control/disease columns
    control_keywords = ['healthy', 'con', 'ref', 'wt', 'norm','c_','ctrl','reduction']
    control_cols = [col for col in expr_cols if any(keyword in col.lower() for keyword in control_keywords)]
    #control_cols = [col for col in expr_cols if col.endswith(('N', 'NSS'))]
    disease_cols = [col for col in expr_cols if col not in control_cols]

    if not control_cols or not disease_cols:
        raise ValueError("Not enough columns for control or disease samples.")

    # Map GeneSymbol to Ensembl_ID, drop rows with missing Ensembl_ID
    expression_df['Ensembl_ID'] = expression_df['GeneSymbol'].map(id_dict.set_index(input_id_type)[output_id_type])
    expression_df.dropna(subset=['Ensembl_ID'], inplace=True)
    ensembl_id_mapped_rows = len(expression_df)

    # Filter numeric columns (control/disease columns)
    numeric_cols = control_cols + disease_cols
    expression_df[numeric_cols] = expression_df[numeric_cols].apply(pd.to_numeric, errors='coerce')

    # Group by Ensembl_ID and average numeric columns (control/disease columns)
    expression_df = expression_df.groupby('Ensembl_ID')[numeric_cols].mean().reset_index()

    # Rename columns for clarity
    control_col_names = [f'control{i+1}' for i in range(len(control_cols))]
    disease_col_names = [f'disease{i+1}' for i in range(len(disease_cols))]
    expression_df.columns = ['Ensembl_ID'] + control_col_names + disease_col_names

    # Save DataFrame to file
    expression_df.to_csv(output_file, index=False)

    # Print summary
    print(f"""
    Control samples: {len(control_cols)}, Disease samples: {len(disease_cols)}
    Initial rows: {initial_rows}, Rows after GeneSymbol filtering: {gene_symbol_filtered_rows}, 
    Rows after Ensembl ID mapping: {ensembl_id_mapped_rows}
    DataFrame saved successfully to {output_file}.
    """)

# Run script if called from command line, otherwise provide arguments manually
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 4:
        dataset_file, id_mapping_file, output_file = sys.argv[1:4]
    else:
        # Use absolute paths for better handling in case of running from different directories
        script_dir = os.path.dirname(os.path.abspath(__file__))
        dataset_file = os.path.join(script_dir, "../../data/cancer/breast_cancer/breast_cancer_GSE9574.txt")
        id_mapping_file = os.path.join(script_dir, "../../data/IDs/SYMBOL_to_ENSEMBL.txt")
        # Example
        output_file = os.path.join(script_dir, "../../data/cancer/breast_cancer/breast_cancer_GSE9574_formatted.csv")

    format_expression_data(dataset_file, id_mapping_file, output_file)

