import pandas as pd

# Define file paths and target pathway
query_file = "../../data/single_cell/CML_E-GEOD-76312_membership.csv"
pathway_file = "../../data/pathways/KEGG/KEGG_2022_membership.tsv" 
target = "hsa05220"

# Load DataFrames
query_df = pd.read_csv(query_file, sep='\t')
pathway_df = pd.read_csv(pathway_file, sep='\t')

# Print columns of query_df to debug
print("Query DataFrame columns:", query_df.columns.tolist())
print("First few rows of query_df:\n", query_df.head())

# Reset the index to include Ensembl_ID as a column
query_df.reset_index(inplace=True)

# Filter for the selected pathway
pathway_filtered = pathway_df[pathway_df["Pathway_Name"] == target]

# Check if pathway_filtered is empty
if pathway_filtered.empty:
    print(f"No data found for pathway: {target}")
else:
    # Merge the DataFrames on 'Ensembl_ID'
    comparison_df = pd.merge(
        pathway_filtered[['Ensembl_ID', "Frequency", 'Overlap_Membership']],
        query_df[['Ensembl_ID', 'Crisp_Membership','Fuzzy_Membership']],
        on='Ensembl_ID',
        how='left'  # Left join keeps all genes from pathway_filtered
    )

    # Display the resulting comparison DataFrame
    print(comparison_df)



pathway_ids = [
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