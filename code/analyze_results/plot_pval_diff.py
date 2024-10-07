import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def load_ora_results(file_path):
    """Load ORA results from a CSV file."""
    return pd.read_csv(file_path)

def plot_p_value_comparison(crisp_file, fuzzy_file, output_dir):
    """Plot p-values from crisp and fuzzy files against each other, coloring points by p-value deviation."""
    
    # Load the results
    crisp_df = load_ora_results(crisp_file)
    fuzzy_df = load_ora_results(fuzzy_file)
    
    # Merge based on 'Pathway_Name'
    merged_df = pd.merge(crisp_df[['Pathway_Name', 'p_value']], 
                         fuzzy_df[['Pathway_Name', 'p_value']], 
                         on='Pathway_Name', suffixes=('_Crisp', '_Fuzzy'))
    
    # Calculate the difference in p-values
    merged_df['p_value_diff'] = merged_df['p_value_Crisp'] - merged_df['p_value_Fuzzy']
    
    # Normalize the p-value differences for coloring
    norm = plt.Normalize(merged_df['p_value_diff'].min(), merged_df['p_value_diff'].max())
    colors = plt.cm.viridis(norm(merged_df['p_value_diff']))  # Choose colormap
    
    # Scatter plot: Crisp P-values vs Fuzzy P-values
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(merged_df['p_value_Crisp'], merged_df['p_value_Fuzzy'], 
                          alpha=0.6, color=colors)
    
    # Color bar for reference
    cbar = plt.colorbar(scatter)
    cbar.set_label('P-Value Deviation')
    
    # Plot settings
    plt.title('Crisp vs Fuzzy P-Value Comparison')
    plt.xlabel('Crisp P-Value')
    plt.ylabel('Fuzzy P-Value')
    
    # Set axis range limits
    plt.xlim(1e-15, 1.05)  # Custom range for x-axis (Crisp p-values)
    plt.ylim(1e-15, 1.05)  # Custom range for y-axis (Fuzzy p-values)
    
    # Reference line for equal p-values
    plt.plot([1e-15, 1.05], [1e-15, 1.05], 'r--')

    # Save the scatter plot
    scatter_plot_path = os.path.join(output_dir, 'scatter_plot.png')
    plt.grid(True)
    plt.savefig(scatter_plot_path)
    plt.close()  # Close the figure to free up memory

    # Bar chart of the p-value differences
    plt.figure(figsize=(8, 6))
    merged_df_sorted = merged_df.sort_values('p_value_diff')
    bars = plt.bar(merged_df_sorted['Pathway_Name'], merged_df_sorted['p_value_diff'], color='blue')
    
    # Add value labels on top of each bar
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, yval, round(yval, 10), 
                 ha='center', va='bottom', rotation=90, fontsize=10)
    
    plt.title('Pathways by P-Value Difference')
    plt.xlabel('Pathway Name')
    plt.ylabel('P-Value Difference')
    plt.xticks(rotation=45, ha='right')  # Rotate x-axis labels for better visibility
    plt.grid(axis='y')
    plt.tight_layout()  # Adjust layout to make room for x-axis labels
    
    # Save the bar chart
    bar_chart_path = os.path.join(output_dir, 'bar_chart.png')
    plt.savefig(bar_chart_path)
    plt.close()  # Close the figure to free up memory

# Define file paths
crisp_file = '/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/ORA_output/infection/HIV/HIV_E-GEOD-111727/standard_ora_results_HIV_E-GEOD-111727_Crisp_Membership_Crisp_Membership.csv'
fuzzy_file = '/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/ORA_output/infection/HIV/HIV_E-GEOD-111727/ora_results_HIV_E-GEOD-111727_Crisp_Membership_Overlap_Membership.csv' 

# Define output directory
output_dir = '/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/ORA_output/infection/HIV/HIV_E-GEOD-111727'

# Plot the p-values
plot_p_value_comparison(crisp_file, fuzzy_file, output_dir)
