import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

# Define pathway groups
cancer_pathway_ids = [
    "hsa05210",  # Colorectal cancer
    "hsa05212",  # Pancreatic cancer
    "hsa05225",  # Hepatocellular carcinoma
    "hsa05226",  # Gastric cancer
    "hsa05214",  # Glioma
    "hsa05216",  # Thyroid cancer
    "hsa05221",  # Acute myeloid leukemia
    "hsa05220",  # Chronic myeloid leukemia
    "hsa05217",  # Basal cell carcinoma
    "hsa05218",  # Melanoma
    "hsa05211",  # Renal cell carcinoma
    "hsa05219",  # Bladder cancer
    "hsa05215",  # Prostate cancer
    "hsa05213",  # Endometrial cancer
    "hsa05224",  # Breast cancer
    "hsa05222",  # Small cell lung cancer
    "hsa05223"   # Non-small cell lung cancer
]

viral_pathway_ids = [
    "hsa05166",  # Human T-cell leukemia virus 1 infection
    "hsa05170",  # Human immunodeficiency virus 1 infection
    "hsa05161",  # Hepatitis B
    "hsa05160",  # Hepatitis C
    "hsa05171",  # Coronavirus disease - COVID-19
    "hsa05164",  # Influenza A
    "hsa05162",  # Measles
    "hsa05168",  # Herpes simplex virus 1 infection
    "hsa05163",  # Human cytomegalovirus infection
    "hsa05167",  # Kaposi sarcoma-associated herpesvirus infection
    "hsa05169",  # Epstein-Barr virus infection
    "hsa05165"   # Human papillomavirus infection
]

infection_pathways = [
    'hsa05166',  # Human T-cell leukemia virus 1 infection
    'hsa05170',  # Human immunodeficiency virus 1 infection
    'hsa05161',  # Hepatitis B
    'hsa05160',  # Hepatitis C
    'hsa05171',  # Coronavirus disease - COVID-19
    'hsa05164',  # Influenza A
    'hsa05162',  # Measles
    'hsa05168',  # Herpes simplex virus 1 infection
    'hsa05163',  # Human cytomegalovirus infection
    'hsa05167',  # Kaposi sarcoma-associated herpesvirus infection
    'hsa05169',  # Epstein-Barr virus infection
    'hsa05165',  # Human papillomavirus infection
    'hsa05110',  # Vibrio cholerae infection
    'hsa05120',  # Epithelial cell signaling in Helicobacter pylori infection
    'hsa05130',  # Pathogenic Escherichia coli infection
    'hsa05132',  # Salmonella infection
    'hsa05131',  # Shigellosis
    'hsa05135',  # Yersinia infection
    'hsa05133',  # Pertussis
    'hsa05134',  # Legionellosis
    'hsa05150',  # Staphylococcus aureus infection
    'hsa05152',  # Tuberculosis
    'hsa05146',  # Amoebiasis
    'hsa05144',  # Malaria
    'hsa05145',  # Toxoplasmosis
    'hsa05140',  # Leishmaniasis
    'hsa05142',  # Chagas disease
    'hsa05143'   # African trypanosomiasis
]


# Define file paths
query_file = "../../data/single_cell/HIV_E-GEOD-111727_membership.csv"
pathway_file = "../../data/pathways/KEGG/KEGG_2022_membership.tsv"
base_output_folder = "../../data/ORA_Output/infection/HIV/HIV_E-GEOD-111727/pathway_histograms"

# Load DataFrames
query_df = pd.read_csv(query_file, sep='\t')
pathway_df = pd.read_csv(pathway_file, sep='\t')

# Membership types to loop over in the pathway DataFrame
membership_types = ['Overlap_Membership_max', 'Overlap_Membership','Overlap_Membership_linear']  # Add other membership types if needed

# Loop through each target pathway
for target in infection_pathways:
    # Filter for the selected pathway
    pathway_filtered = pathway_df[pathway_df["Pathway_Name"] == target]

    if pathway_filtered.empty:
        print(f"No data found for pathway: {target}")
        continue

    # Loop over each membership type in pathway_df
    for membership_type in membership_types:
        # Skip if the membership column is not present in the DataFrame
        if membership_type not in pathway_filtered.columns:
            print(f"Membership type {membership_type} not found in pathway file for {target}")
            continue

        # Create output directory for this membership type if it doesn't exist
        output_folder = os.path.join(base_output_folder, membership_type)
        os.makedirs(output_folder, exist_ok=True)

        # Merge the DataFrames on 'Ensembl_ID'
        comparison_df = pd.merge(
            pathway_filtered[['Ensembl_ID', "Frequency", membership_type]],
            query_df[['Ensembl_ID', 'Crisp_Membership']],
            on='Ensembl_ID',
            how='left'
        )

        # Filter the overlap memberships for Crisp Membership 0 and 1
        overlap_memberships_crisp_0 = comparison_df[comparison_df['Crisp_Membership'] == 0][membership_type].dropna()
        overlap_memberships_crisp_1 = comparison_df[comparison_df['Crisp_Membership'] == 1][membership_type].dropna()

        # Create a histogram with stacking
        plt.figure(figsize=(10, 6))
        bins = np.linspace(0, 1, 31)  # 30 bins from 0 to 1
        
        # Create the second histogram for Crisp Membership 1, stacking on top of the first
        hist_0, bins_0 = np.histogram(overlap_memberships_crisp_1, bins=bins)
        plt.bar(bins_0[:-1], hist_0, width=np.diff(bins), color='red', edgecolor='black', alpha=0.7, label='Crisp Membership = 1')

        # Create the first histogram for Crisp Membership 0
        hist_1, bins_1 = np.histogram(overlap_memberships_crisp_0, bins=bins)
        plt.bar(bins_1[:-1], hist_1, width=np.diff(bins), bottom=hist_0, color='blue', edgecolor='black', alpha=0.7, label='Crisp Membership = 0')

        # Customize the plot
        pathway_description = pathway_filtered['Description'].values[0]  # Assuming there's a column for description
        plt.title(f'Stacked Histogram of {membership_type} for {target} - {pathway_description}')
        plt.xlabel(f'{membership_type}')
        plt.ylabel('Frequency')
        plt.grid(axis='y', alpha=0.75)
        plt.xticks(np.arange(0, 1.1, step=0.1))  # Set x-ticks based on the data range
        plt.legend()  # Show the legend

        # Save the plot to the membership-specific folder
        output_file = os.path.join(output_folder, f"{target}_{pathway_description}_{membership_type}_histogram.png")
        plt.savefig(output_file)  # Save the figure
        plt.close()  # Close the plot to avoid displaying it in an interactive session

        print(f"Stacked plot saved for {membership_type} in pathway {target} to {output_file}")

