import pandas as pd
import matplotlib.pyplot as plt
import os

def analyze_expansion(expanded_pathways_path, output_dir):
    try:
        # Read the expanded pathways file
        pathways_expanded = pd.read_csv(expanded_pathways_path, sep="\t")
        
        # Ensure the output directory exists
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Define normalization function
        def normalize(value):
            min_value = 0.4
            max_value = 1
            if value < min_value:
                return 0
            elif value > max_value:
                return 1
            else:
                return (value - min_value) / (max_value - min_value)
        
        # Group pathways by 'Pathway_Name' column
        grouped_pathways = pathways_expanded.groupby('Pathway_Name')
        
        # Create histograms for each pathway
        for pathway_name, group_df in grouped_pathways:
            # Ensure the figure is closed before creating new ones
            plt.close('all')
            
            # Create and save histogram for original Expansion_Membership values
            plt.figure(figsize=(10, 6))
            plt.hist(group_df['Expansion_Membership'].dropna(), bins=30, edgecolor='black', density=False)
            plt.title(f'Histogram of Expansion_Membership for {pathway_name}')
            plt.xlabel('Expansion_Membership')
            plt.ylabel('Frequency')
            plt.savefig(os.path.join(output_dir, f'{pathway_name}_original_histogram.png'))
            plt.close()
            
            # Filter and create histogram for Expansion_Membership values above 0.4
            filtered_df = group_df[group_df['Expansion_Membership'] > 0.4]
            plt.figure(figsize=(10, 6))
            plt.hist(filtered_df['Expansion_Membership'].dropna(), bins=30, edgecolor='black', density=False)
            plt.title(f'Histogram of Expansion_Membership (filtered) for {pathway_name}')
            plt.xlabel('Expansion_Membership')
            plt.ylabel('Frequency')
            plt.savefig(os.path.join(output_dir, f'{pathway_name}_filtered_histogram.png'))
            plt.close()
            
            # Normalize and create histogram for normalized Expansion_Membership values
            filtered_df['Normalized_Membership'] = filtered_df['Expansion_Membership'].apply(normalize)
            plt.figure(figsize=(10, 6))
            plt.hist(filtered_df['Normalized_Membership'].dropna(), bins=30, edgecolor='black', density=False)
            plt.title(f'Histogram of Normalized Expansion_Membership (filtered) for {pathway_name}')
            plt.xlabel('Normalized_Expansion_Membership')
            plt.ylabel('Frequency')
            plt.savefig(os.path.join(output_dir, f'{pathway_name}_normalized_histogram.png'))
            plt.close()
            
            print(f"Histograms for {pathway_name} saved to {output_dir}")
    
    except pd.errors.ParserError as e:
        print("Parsing Error:", e)
    except Exception as e:
        print("An error occurred:", e)

# Example usage
expanded_pathways_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/KEGG_2022_entrez_with_membership.tsv"
output_dir = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/Expansion_Histograms"

analyze_expansion(expanded_pathways_path, output_dir)
