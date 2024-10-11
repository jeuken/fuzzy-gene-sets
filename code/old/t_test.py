import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

def calculate_statistics_from_expression_dataframe(expression_path, output_path):
    """
    Calculates T-statistics and p-values for each gene based on expression data,
    ranks genes by T-statistic from high to low, adjusts p-values, and saves the results to a specified output file.
    
    Args:
        expression_path (str): Path to the expression data file.
        output_path (str): Path where the results will be saved.
    """
    try:
        # Load the expression DataFrame
        df = pd.read_csv(expression_path, sep=',')  # Adjust separator as needed
        
        # Extract column names corresponding to control and disease conditions
        control_data = [col for col in df.columns if col.startswith('control')]
        disease_cols = [col for col in df.columns if col.startswith('disease')]
        
        # Convert expression data columns to numeric with float dtype, forcing errors to NaN
        df[control_data + disease_cols] = df[control_data + disease_cols].apply(pd.to_numeric, errors='coerce')
        
        # Drop rows where 'Ensembl_ID' is missing
        initial_row_count = len(df)
        df.dropna(subset=['Ensembl_ID'], inplace=True)
        missing_id_count = initial_row_count - len(df)
        print(f"Rows dropped due to missing Ensembl IDs: {missing_id_count}")
        
        # Count NaNs before T-test
        total_nans_control = df[control_data].isna().sum().sum()
        total_nans_disease = df[disease_cols].isna().sum().sum()
        total_nans = total_nans_control + total_nans_disease
        print(f"Total NaNs to be ignored: {total_nans}")
        
        # Calculate T-statistics and p-values using vectorized operations
        control_data_values = df[control_data].values
        disease_data_values = df[disease_cols].values
        
        # Calculate T-statistics and p-values for each gene, omitting NaN values for each gene
        t_stats, p_vals = ttest_ind(control_data_values, disease_data_values, axis=1, equal_var=False, nan_policy='omit')
        
        # Adjust p-values using Benjamini-Hochberg procedure
        adjusted_p_vals = multipletests(p_vals, method='fdr_bh')[1]
        
        # Create a DataFrame with the results and correct column names
        results_df = pd.DataFrame({
            'Ensembl_ID': df['Ensembl_ID'],
            't': t_stats,
            'p_value': p_vals,
            'adj_p_value': adjusted_p_vals
        })
        
        # Rank genes by T-statistic from high to low
        results_df = results_df.sort_values(by='t', ascending=False)
        
        # Save the results DataFrame to a file
        results_df.to_csv(output_path, sep=',', index=False)  # Adjust separator as needed
        
        print(f"Results saved successfully to {output_path}.")
        
        # Display the first few rows of the results DataFrame
        print(results_df.head())
        
    except Exception as e:
        print(f"An error occurred: {e}")

# Set the input and output paths
input_path = "../../data/gemma/Alzheimers_GSE95587_formatted.csv"
output_path = "../../data/gemma/Alzheimers_GSE95587_t_stats.csv"

# Execute the function
calculate_statistics_from_expression_dataframe(input_path, output_path)

