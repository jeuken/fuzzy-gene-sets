import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

def calculate_statistics(expression_path, output_path):
    # Load the expression DataFrame
    df = pd.read_csv(expression_path, sep=',')
    
    # Identify control and disease columns
    control_cols = [col for col in df.columns if col.startswith('control')]
    disease_cols = [col for col in df.columns if col.startswith('disease')]
    
    # Drop rows with missing Ensembl_ID and convert expression data to numeric
    df.dropna(subset=['Ensembl_ID'], inplace=True)
    df[control_cols + disease_cols] = df[control_cols + disease_cols].apply(pd.to_numeric, errors='coerce')

    # Perform T-tests
    t_stats, p_vals = ttest_ind(df[control_cols].values, df[disease_cols].values, axis=1, equal_var=False, nan_policy='omit')
    
    # Adjust p-values using Benjamini-Hochberg correction
    adjusted_p_vals = multipletests(p_vals, method='fdr_bh')[1]
    
    # Create and rank results by T-statistic
    results_df = pd.DataFrame({'Ensembl_ID': df['Ensembl_ID'], 't': t_stats, 'p_value': p_vals, 'adj_p_value': adjusted_p_vals})
    results_df.sort_values(by='t', ascending=False, inplace=True)
    
    # Count significant genes (adj_p_value <= 0.05)
    significant_genes = (results_df['adj_p_value'] <= 0.05).sum()
    total_genes = len(results_df)
    percentage_significant = (significant_genes / total_genes) * 100
    print(f"Number of significant genes (adj_p_value <= 0.05): {significant_genes}")
    print(f"Percentage of significant genes: {percentage_significant:.2f}%")
    
    # Save the results
    results_df.to_csv(output_path, sep=',', index=False)
    print(f"Results saved to {output_path}.")
    return results_df

# Command-line or direct execution
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 3:
        input_path, output_path = sys.argv[1:3]
    else:
        input_path = "../../data/cancer/breast_cancer/breast_cancer_GSE9574_formatted.csv"
        output_path = "../../data/cancer/breast_cancer/breast_cancer_GSE9574_t_stats.csv"
    
t_stats = calculate_statistics(input_path, output_path)
