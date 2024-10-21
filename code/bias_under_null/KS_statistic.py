import pandas as pd
from scipy import stats

def load_p_values(file_path):
    """Load p-values from a CSV file."""
    return pd.read_csv(file_path, sep='\t')

def ks_test_all_p_values(p_values):
    """Perform K-S test on all p-values."""
    ks_statistic, p_value = stats.kstest(p_values, 'uniform')
    return ks_statistic, p_value

def ks_test_per_pathway(p_values_df):
    """Perform K-S test for each pathway."""
    results = {}
    for pathway in p_values_df.index:
        pathway_p_values = p_values_df.loc[pathway].dropna().values
        if len(pathway_p_values) > 0:
            ks_statistic, p_value = stats.kstest(pathway_p_values, 'uniform')
            results[pathway] = {'KS_Statistic': ks_statistic, 'p_value': p_value}
        else:
            results[pathway] = {'KS_Statistic': None, 'p_value': None}
    return results

def determine_uniformity(p_values_file, output_file):
    # Load p-values
    p_values_df = load_p_values(p_values_file)

    # Convert all values to numeric, forcing errors to NaN
    p_values_df = p_values_df.apply(pd.to_numeric, errors='coerce')

    # K-S test on all p-values
    all_p_values = p_values_df.values.flatten()
    all_p_values = all_p_values[~pd.isnull(all_p_values)]  # Remove NaN values
    ks_stat_all, p_value_all = ks_test_all_p_values(all_p_values)

    # K-S test per pathway
    ks_results_per_pathway = ks_test_per_pathway(p_values_df)

    # Create summary DataFrame
    summary_df = pd.DataFrame({
        'Pathway_Name': ks_results_per_pathway.keys(),
        'KS_Statistic': [result['KS_Statistic'] for result in ks_results_per_pathway.values()],
        'p_value': [result['p_value'] for result in ks_results_per_pathway.values()]
    })

    # Add overall K-S test results
    summary_df.loc[len(summary_df)] = ['Overall', ks_stat_all, p_value_all]

    # Save results to CSV
    summary_df.to_csv(output_file, index=False)

    print(f"K-S test results saved to: {output_file}")

# Define paths
p_values_file = '../../data/bias_under_null/null_distributions_crisp_overlap.csv'  # Input file with p-values
output_file = '../../data/bias_under_null/ks_test_results_crisp_overlap.csv'  # Output file for K-S test results

# Run the main function
determine_uniformity(p_values_file, output_file)
