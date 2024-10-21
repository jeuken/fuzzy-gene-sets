import pandas as pd
import numpy as np
from scipy import stats
import os
import matplotlib.pyplot as plt

def load_p_values(file_path):
    """Load p-values from a CSV file."""
    return pd.read_csv(file_path, sep='\t')

def chi_square_test(p_values, bins=10):
    """Perform Chi-square goodness of fit test for uniform distribution."""
    # Define the bin edges for uniform distribution
    observed_freq, bin_edges = np.histogram(p_values, bins=bins, range=(0, 1))

    # Expected frequencies assuming uniform distribution
    expected_freq = np.ones_like(observed_freq) * len(p_values) / bins

    # Perform the Chi-square test
    chi_statistic, p_value = stats.chisquare(f_obs=observed_freq, f_exp=expected_freq)
    
    return chi_statistic, p_value, observed_freq, expected_freq, bin_edges

def chi_square_test_per_pathway(p_values_df, bins=10):
    """Perform Chi-square test for each pathway."""
    results = {}
    for pathway in p_values_df.index:
        pathway_p_values = p_values_df.loc[pathway].dropna().values
        if len(pathway_p_values) > 0:
            chi_stat, p_value, observed, expected, bin_edges = chi_square_test(pathway_p_values, bins)
            results[pathway] = {
                'Chi_Statistic': chi_stat,
                'p_value': p_value,
                'Observed_Frequencies': observed,
                'Expected_Frequencies': expected,
                'Bin_Edges': bin_edges
            }
        else:
            results[pathway] = {'Chi_Statistic': None, 'p_value': None}
    return results

def plot_histogram_with_fit(pathway_name, observed, expected, bin_edges, output_dir):
    """Plot histogram with the observed and expected frequencies for each pathway."""
    plt.figure(figsize=(10, 6))
    plt.hist(bin_edges[:-1], bin_edges, weights=observed, alpha=0.7, color='blue', label='Observed')
    plt.plot(bin_edges[:-1] + (bin_edges[1] - bin_edges[0]) / 2, expected, 'r--', label='Expected (Uniform)')
    
    plt.title(f'Chi-Square Test for {pathway_name}')
    plt.xlabel('p-value')
    plt.ylabel('Frequency')
    plt.legend()
    plt.grid(True)
    
    # Save the plot
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, f'{pathway_name}_chi_square.png'))
    plt.close()

def determine_uniformity(p_values_file, output_file, bins=10):
    # Load p-values
    p_values_df = load_p_values(p_values_file)

    # Convert all values to numeric, forcing errors to NaN
    p_values_df = p_values_df.apply(pd.to_numeric, errors='coerce')

    # Perform Chi-square test on all p-values
    all_p_values = p_values_df.values.flatten()
    all_p_values = all_p_values[~pd.isnull(all_p_values)]  # Remove NaN values
    chi_stat_all, p_value_all, observed_all, expected_all, bin_edges_all = chi_square_test(all_p_values, bins)

    # Perform Chi-square test per pathway
    chi_results_per_pathway = chi_square_test_per_pathway(p_values_df, bins)

    # Create summary DataFrame
    summary_df = pd.DataFrame({
        'Pathway_Name': chi_results_per_pathway.keys(),
        'Chi_Statistic': [result['Chi_Statistic'] for result in chi_results_per_pathway.values()],
        'p_value': [result['p_value'] for result in chi_results_per_pathway.values()]
    })

    # Add overall Chi-square test results
    summary_df.loc[len(summary_df)] = ['Overall', chi_stat_all, p_value_all]

    # Save results to CSV
    summary_df.to_csv(output_file, index=False)

    print(f"Chi-Square test results saved to: {output_file}")

    # Save histograms with observed vs expected frequencies for all pathways
    output_dir = '../../data/bias_under_null/chi_square_plots/'  # Directory to save the plots
    plot_histogram_with_fit('Overall', observed_all, expected_all, bin_edges_all, output_dir)

    for pathway, result in chi_results_per_pathway.items():
        if result['Chi_Statistic'] is not None:
            plot_histogram_with_fit(pathway, result['Observed_Frequencies'], result['Expected_Frequencies'], result['Bin_Edges'], output_dir)

# Define paths
p_values_file = '../../data/bias_under_null/null_distributions.csv'  # Input file with p-values
output_file = '../../data/bias_under_null/chi_square_test_results.csv'  # Output file for Chi-square test results

# Run the main function
determine_uniformity(p_values_file, output_file, bins=10)

