import pandas as pd
import os

def filter_and_normalize(expanded_pathways_path, output_file):
    try:
        # Read the expanded pathways file
        pathways_expanded = pd.read_csv(expanded_pathways_path, sep="\t")
        
        # Filter out rows with Expansion_Membership below 0.4
        pathways_expanded = pathways_expanded[pathways_expanded['Expansion_Membership'] >= 0.4]
        
        # Define normalization function
        def normalize(value):
            min_value = 0.4
            max_value = 1
            return (value - min_value) / (max_value - min_value)
        
        # Apply normalization to the remaining rows
        pathways_expanded['Expansion_Membership'] = pathways_expanded['Expansion_Membership'].apply(lambda x: normalize(x) if pd.notna(x) else x)
      
        # Save the new DataFrame to a file
        pathways_expanded.to_csv(output_file, sep="\t", index=False)
        
        print(f"Filtered and normalized DataFrame saved to {output_file}")
    
    except pd.errors.ParserError as e:
        print("Parsing Error:", e)
    except Exception as e:
        print("An error occurred:", e)

# Example usage
expanded_pathways_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/expanded_mapped_pathways.tsv"
output_file = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/pathways_with_membership.tsv"

filter_and_normalize(expanded_pathways_path, output_file)
