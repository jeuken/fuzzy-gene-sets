import pandas as pd  # Import pandas for data manipulation and analysis
import matplotlib.pyplot as plt  # Import matplotlib for plotting
import os  # Import os for file path operations
import numpy as np  # Import numpy for numerical operations

# Step 1: Load the data
null_distribution_path = '../../../data/bias_under_null/null_distributions_crisp_overlap_dp.csv'  # Path to the null distribution data
results = pd.read_csv(null_distribution_path, sep='\t', index_col=0)  # Load the data, using the first column as the index

# Ensure the output directory exists before saving
output_dir_hist = '../../../data/bias_under_null/null_distributions/histograms_crisp_overlap_dp/'  # Path to save the histograms
os.makedirs(output_dir_hist, exist_ok=True)  # Create the directory if it does not exist

# Ensure the output directory exists before saving
output_dir_cdf = '../../../data/bias_under_null/null_distributions/cdf_crisp_overlap_dp/'  # Path to save the histograms
os.makedirs(output_dir_hist, exist_ok=True)  # C

# Ensure the output directories exist before saving
output_dir_pathway_hist = '../../../data/bias_under_null/null_distributions/histograms_crisp_overlap_dp/'  # Path to save the histograms
output_dir_pathway_cdf = '../../../data/bias_under_null/null_distributions/cdf_crisp_overlap_dp/'  # Path to save the CDFs
os.makedirs(output_dir_pathway_hist, exist_ok=True) 
os.makedirs(output_dir_pathway_cdf, exist_ok=True) 


# Step 2: Reshape the DataFrame for easier access to p-values
# Extract column names that contain '_p_value' for further analysis
p_value_columns = [col for col in results.columns if '_p_value' in col]

# Step 3: Create Overall Histogram
plt.figure(figsize=(10, 6))
plt.hist(results[p_value_columns].values.flatten(), bins=100, alpha=0.7, color='blue')
plt.title('Overall Histogram of p-values')  # Set the title of the histogram
plt.xlabel('p-value')  # Label for the x-axis
plt.ylabel('Frequency')  # Label for the y-axis
plt.grid(axis='y')  # Add grid lines along the y-axis for better readability
plt.tight_layout()  # Automatically adjust subplot parameters for a clean layout


# Save the overall histogram
plt.savefig(os.path.join(output_dir_hist, 'overall_p_value_histogram.png'))  # Save the histogram to a file
plt.close()  # Close the figure to free up memory

# Step 4: Create Overall CDF
plt.figure(figsize=(10, 6))
# Flatten the array of p-values and calculate ECDF
all_p_values = results[p_value_columns].values.flatten()
sorted_all_p_values = np.sort(all_p_values)
ecdf_all = np.arange(1, len(sorted_all_p_values) + 1) / len(sorted_all_p_values)

plt.step(sorted_all_p_values, ecdf_all, where='post', color='green', label='Empirical CDF')
plt.plot([0, 1], [0, 1], color='red', linestyle='--', label='Expected CDF (Uniform)')
plt.title('Overall ECDF of p-values')  # Set the title
plt.xlabel('p-value')  # Label for the x-axis
plt.ylabel('Cumulative Probability')  # Label for the y-axis
plt.grid()  # Add grid lines
plt.legend()  # Add a legend
plt.tight_layout()  # Automatically adjust subplot parameters for a clean layout

# Save the overall ECDF plot
plt.savefig(os.path.join(output_dir_cdf, 'overall_ecdf_p_values.png'))  # Save the ECDF plot to a file
plt.close()  # Close the figure to free up memory

# Step 5: Create Histograms and CDFs per Pathway
# Use the index for pathway names, which allows for individual pathway access
pathway_names = results.index  # Access pathway names from the DataFrame index


# Iterate through each pathway to create individual histograms and ECDFs
for pathway in pathway_names:
    # Extract p-values for the current pathway from the DataFrame
    pathway_p_values = results.loc[pathway, p_value_columns].values.flatten()
    
    # Step 6: Plot Histogram for the current pathway
    plt.figure(figsize=(10, 6))
    plt.hist(pathway_p_values, bins=100, alpha=0.7, color='blue')
    plt.title(f'Histogram of p-values for {pathway}')  # Set the title with the pathway name
    plt.xlabel('p-value')  # Label for the x-axis
    plt.ylabel('Frequency')  # Label for the y-axis
    plt.grid(axis='y')  # Add grid lines along the y-axis
    plt.tight_layout()  # Automatically adjust subplot parameters for a clean layout

    # Save the histogram for the pathway
    plt.savefig(os.path.join(output_dir_pathway_hist, f'{pathway}_histogram.png'))  # Save the histogram to a file
    plt.close()  # Close the figure to free up memory

    # Step 7: Plot ECDF for the current pathway
    plt.figure(figsize=(10, 6))
    
    # Sort the p-values for ECDF
    sorted_pathway_p_values = np.sort(pathway_p_values)
    ecdf_pathway = np.arange(1, len(sorted_pathway_p_values) + 1) / len(sorted_pathway_p_values)  # Compute ECDF
    
    # Plot the ECDF
    plt.step(sorted_pathway_p_values, ecdf_pathway, where='post', color='green', label='Empirical CDF')
    
    # Add the expected CDF reference line
    plt.plot([0, 1], [0, 1], color='red', linestyle='--', label='Expected CDF (Uniform)')

    plt.title(f'ECDF of p-values for {pathway}')  # Set the title with the pathway name
    plt.xlabel('p-value')  # Label for the x-axis
    plt.ylabel('Cumulative Probability')  # Label for the y-axis
    plt.grid()  # Add grid lines
    plt.legend()  # Add a legend
    plt.tight_layout()  # Automatically adjust subplot parameters for a clean layout

    # Save the ECDF plot for each pathway
    plt.savefig(os.path.join(output_dir_pathway_cdf, f'{pathway}_ecdf_p_values.png'))  # Save the ECDF plot to a file
    plt.close()  # Close the figure to avoid display issues and free up memory
