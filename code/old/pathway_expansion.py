import pandas as pd
import numpy as np

def pathway_expansion(string_scores_path, pathways_path, expanded_pathways_path):
    try:
        # Read the CSV files with tab-separated values
        string = pd.read_csv(string_scores_path, sep="\t")
        pathways = pd.read_csv(pathways_path, sep="\t")
        
        # Add a new column "Expansion_Membership" and initialize with 1s for existing genes
        pathways["Expansion_Membership"] = 1
        
        # Group pathways by 'Pathway_Name' column
        grouped_pathways = pathways.groupby('Pathway_Name')
        
        # Create a list to store new rows that need to be added
        new_rows = []
        
        # Initialize a counter for completed pathways
        completed_pathways_count = 0

        # Iterate over each group
        for pathway_name, pathway_df in grouped_pathways:
            pathway_genes = set(pathway_df['Gene_ID'])  # Assuming 'Gene_ID' is the correct column name
            
            # Extract pathway description from the first row (assuming all rows have the same description for the pathway)
            pathway_description = pathway_df['Description'].iloc[0]  # Adjust if 'Description' is the correct column name
            
            # Filter string DataFrame for rows where 'Gene2' is in the pathway genes
            pathway_string = string[string['Gene2'].isin(pathway_genes)]
            
            # Calculate average STRING score for 'Gene1' genes across genes in the pathway
            gene1_scores = pathway_string.groupby('Gene1')['string_score'].mean()  # Assuming 'string_score' is the correct column
            
            # Check if genes in 'Gene1' are not in the pathway and prepare them to be added
            for gene, avg_score in gene1_scores.items():
                if gene not in pathway_genes:
                    # Create a new row dictionary with the same structure as pathways_df
                    new_row = {
                        'Pathway_Name': pathway_name,  # Same pathway name
                        'Gene_ID': gene,               # New gene from Gene1
                        'Expansion_Membership': avg_score / 1000,  # Set Expansion_Membership to avg score divided by 1000
                        'Description': pathway_description  # Copy the pathway description
                    }
                    
                    # Add other columns from the pathways_df, initializing them with NaN if they aren't in the new_row
                    for col in pathways.columns:
                        if col not in new_row:
                            new_row[col] = np.nan  # Set missing columns to NaN
                    
                    # Append the new row to the list of new rows
                    new_rows.append(new_row)
            
            # Increment the completed pathways counter
            completed_pathways_count += 1
            print(f"Number of pathways completed: {completed_pathways_count}")
        
        # Convert new rows to a DataFrame and append to the original pathways DataFrame
        new_rows_df = pd.DataFrame(new_rows)
        pathways_expanded = pd.concat([pathways, new_rows_df], ignore_index=True)
        
        # Save the expanded DataFrame to the specified path
        pathways_expanded.to_csv(expanded_pathways_path, sep="\t", index=False)
        print(f"Expanded pathways saved to {expanded_pathways_path}")
        
        # Calculate and print max and min values below 1 in the Expansion_Membership column
        filtered_df = pathways_expanded[pathways_expanded['Expansion_Membership'] < 1]
        max_value = filtered_df['Expansion_Membership'].max()
        min_value = filtered_df['Expansion_Membership'].min()
        
        print(f"Maximum value in 'Expansion_Membership' (below 1): {max_value}")
        print(f"Minimum value in 'Expansion_Membership' (below 1): {min_value}")
        
        return pathways_expanded

    except pd.errors.ParserError as e:
        print("Parsing Error:", e)
        return None
    except Exception as e:
        print("An error occurred:", e)
        return None

# Example usage
string_scores_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/STRING/STRING_GGI.txt"
pathways_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/mapped_pathways.tsv"
expanded_pathways_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/expanded_mapped_pathways.tsv"

results = pathway_expansion(string_scores_path, pathways_path, expanded_pathways_path)
print(results)

