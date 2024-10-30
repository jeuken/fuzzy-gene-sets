import pandas as pd
import matplotlib.pyplot as plt
import os

def generate_membership_histograms(data, pathway_ids, category):
    # Set up output directory based on category
    base_output_dir = f'/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/data_analysis/overlap_plots/{category}'
    
    for target_pathway in pathway_ids:
        # Create output directory for the current target pathway
        output_dir = os.path.join(base_output_dir, target_pathway)
        os.makedirs(output_dir, exist_ok=True)
        
        # Get genes for the current target pathway
        target_genes = data.loc[data['Pathway_Name'] == target_pathway, 'Ensembl_ID'].unique()
        
        for neg_pathway in pathway_ids:
            if neg_pathway == target_pathway:
                continue  # Skip if the pathway is the same as the target
            
            # Filter data for the negative pathway
            pathway_data = data.loc[data['Pathway_Name'] == neg_pathway]
            
            if pathway_data.empty:
                continue
            
            # Determine genes in the target pathway
            pathway_data.loc[:, 'In_Target'] = pathway_data['Ensembl_ID'].isin(target_genes)
            
            # Separate membership values
            in_target_memberships = pathway_data.loc[pathway_data['In_Target'], 'Overlap_Membership']
            not_in_target_memberships = pathway_data.loc[~pathway_data['In_Target'], 'Overlap_Membership']
            
            # Plot setup
            fig, ax = plt.subplots(figsize=(10, 6))
            
            # Create the histogram
            ax.hist(
                [in_target_memberships, not_in_target_memberships],
                bins=30,
                stacked=True,
                color=['red', 'blue'],
                label=['In Target Pathway', 'Not in Target Pathway']
            )
            
            # Add title with pathway descriptions
            neg_description = pathway_data['Description'].iloc[0]
            target_description = data.loc[data['Pathway_Name'] == target_pathway, 'Description'].iloc[0]
            ax.set_title(f"{neg_pathway} - {neg_description} vs {target_pathway} - {target_description}")
            
            # Labels and legend
            ax.set_xlabel('Membership Value')
            ax.set_ylabel('Gene Count')
            ax.legend()
            
            # Save the plot
            plot_filename = os.path.join(output_dir, f"{neg_pathway}_membership_histogram.png")
            plt.tight_layout()
            plt.savefig(plot_filename)
            plt.close()

# Load data
file_path = '/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/KEGG_2022_pathway_memberships.tsv' 
data = pd.read_csv(file_path, sep='\t')

# Define pathway IDs
cancer_pathway_ids = ["hsa05210", "hsa05212", "hsa05225", "hsa05226", "hsa05214", "hsa05216", "hsa05221",
                      "hsa05220", "hsa05217", "hsa05218", "hsa05211", "hsa05219", "hsa05215", "hsa05213",
                      "hsa05224", "hsa05222", "hsa05223"]

infection_pathway_ids = ['hsa05166', 'hsa05170', 'hsa05161', 'hsa05160', 'hsa05171', 'hsa05164', 'hsa05162',
                         'hsa05168', 'hsa05163', 'hsa05167', 'hsa05169', 'hsa05165', 'hsa05110', 'hsa05120',
                         'hsa05130', 'hsa05132', 'hsa05131', 'hsa05135', 'hsa05133', 'hsa05134', 'hsa05150',
                         'hsa05152', 'hsa05146', 'hsa05144', 'hsa05145', 'hsa05140', 'hsa05142', 'hsa05143']

# Run for cancer pathways
generate_membership_histograms(data, cancer_pathway_ids, 'cancer')

# Run for infection pathways
generate_membership_histograms(data, infection_pathway_ids, 'infection')

