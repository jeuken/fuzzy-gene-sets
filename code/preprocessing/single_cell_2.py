import pandas as pd
import numpy as np
from scipy.stats import f_oneway
from joblib import Parallel, delayed, cpu_count
from qvalue import qvalues
from tqdm import tqdm  
import os

# Specify labels and disease name at the top
negative_label = "normal"           # Set the negative label here
positive_label = "HIV infection"    # Set the positive label here
disease_name = "HIV"                # Specify the short disease name for output files

# Function definitions
def anova_proc_scipy(data_groupby, gene_name):
    f, p = f_oneway(*data_groupby[gene_name].apply(list))
    return (gene_name, p, f)

def anova_genes_scipy(data, label, n_jobs):
    data['label'] = label
    data_gb = data.groupby('label', observed=False)
    
    result = Parallel(n_jobs=n_jobs)(
        delayed(anova_proc_scipy)(data_gb, gene) for gene in tqdm(data.columns[:-1], desc="Running ANOVA", unit="gene")
    )
    
    result_df = pd.DataFrame(data=result, columns=['Ensembl_ID', 'p', 'f']).set_index('Ensembl_ID')
    return result_df

def log_fold_change_proc(data_groupby, gene_name, order, positive_label):
    unpacked = data_groupby[gene_name].apply(list)
    if order[0] == positive_label:
        mean1 = np.mean(unpacked[order[0]])  # Positive label
        mean2 = np.mean(unpacked[order[1]])  # Negative label
    else:
        mean1 = np.mean(unpacked[order[1]])  # Positive label
        mean2 = np.mean(unpacked[order[0]])  # Negative label
    fc = mean1 - mean2
    return (gene_name, fc)

def log_fold_change(data, label, n_jobs, positive_label):
    data['label'] = label
    data_gb = data.groupby('label', observed=False)
    order = label.cat.categories  # Get the order of categories
    
    result = Parallel(n_jobs=n_jobs)(
        delayed(log_fold_change_proc)(data_gb, gene, order, positive_label) for gene in tqdm(data.columns[:-1], desc="Calculating Log Fold Change", unit="gene")
    )
    
    result_df = pd.DataFrame(data=result, columns=['Ensembl_ID', 'logfc']).set_index('Ensembl_ID')
    return result_df

# Main script - Part 2
# Load the saved data from Part 1
metadata = pd.read_pickle("../../data/single_cell/metadata.pkl")
adata = pd.read_pickle("../../data/single_cell/adata.pkl")

# Ensure that label is categorical
metadata['Factor Value[disease]'] = metadata['Factor Value[disease]'].astype('category')

# Filter the data based on the selected labels
labels_of_interest = [negative_label, positive_label]
metadata_filtered = metadata[metadata['Factor Value[disease]'].isin(labels_of_interest)]
adata_filtered = adata.loc[metadata_filtered.index]
label_filtered = metadata_filtered['Factor Value[disease]']

# Run ANOVA across genes
n_jobs = cpu_count() - 1
results = anova_genes_scipy(adata_filtered, label_filtered, n_jobs)

# Apply q-value correction
results = qvalues(results)

# Run log fold change calculation
logfc = log_fold_change(adata_filtered, label_filtered, 10, positive_label)
results['logfc'] = logfc['logfc']

# Define output path and save the results
output_path = "../../data/single_cell"
disease = disease_name  # Use the short disease name for saving results

# Ensure the output directory exists
os.makedirs(output_path, exist_ok=True)

# Save the results DataFrame to a CSV file
study = 'E-GEOD-111727'
results.to_csv(os.path.join(output_path, f'{disease}_{study}.csv'), index=True)

print(f"Results saved to {os.path.join(output_path, f'{disease}_{study}.csv')}")
