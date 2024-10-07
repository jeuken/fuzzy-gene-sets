import scanpy as sc
import os

# Main script - Part 1
study = 'E-GEOD-111727'

# Load dataset from Expression Atlas
adata = sc.datasets.ebi_expression_atlas(study)

# Basic filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Library size correction
sc.pp.normalize_total(adata, target_sum=1e4)

# Logarithmize the data
sc.pp.log1p(adata)

# Get the metadata and label
metadata = adata.obs
adata = adata.to_df()

# Print the available disease labels
possible_labels = metadata['Factor Value[disease]'].unique()
print("Possible disease labels:", possible_labels)

# Save metadata and expression data for use in Part 2
metadata_path = "../../data/single_cell/metadata.pkl"
adata_path = "../../data/single_cell/adata.pkl"

# Ensure the output directory exists
os.makedirs("../../data/single_cell", exist_ok=True)

# Save metadata and processed data for later use
metadata.to_pickle(metadata_path)
adata.to_pickle(adata_path)

print("Part 1 complete. The metadata and adata have been saved.")
