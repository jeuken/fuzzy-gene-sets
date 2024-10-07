import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import f_oneway
from joblib import Parallel, delayed, cpu_count
from qvalue import qvalues
import os
from tqdm import tqdm  # Import tqdm

def anova_proc_scipy(data_groupby, gene_name):
    f, p = f_oneway(*data_groupby[gene_name].apply(list))
    return (gene_name, p, f)

def anova_genes_scipy(data, label, n_jobs):
    data['label'] = label 
    data_gb = data.groupby('label', observed=False)
    
    # Use tqdm to add a progress bar
    result = Parallel(n_jobs=n_jobs)(
        delayed(anova_proc_scipy)(data_gb, gene) for gene in tqdm(data.columns[:-1], desc="Running ANOVA", unit="gene")
    )
    
    result_df = pd.DataFrame(data=result, columns=['Ensembl_ID', 'p', 'f']).set_index('Ensembl_ID')
    return result_df

def log_fold_change_proc(data_groupby, gene_name, order):
    unpacked = data_groupby[gene_name].apply(list)
    mean1 = np.mean(unpacked[order[0]])
    mean2 = np.mean(unpacked[order[1]])
    fc = mean1 - mean2
    return (gene_name, fc)

def log_fold_change(data, label, n_jobs):
    data['label'] = label 
    data_gb = data.groupby('label', observed=False)
    order = data_gb[data.columns[0]].apply(list).index
    result = Parallel(n_jobs=n_jobs)(
        delayed(log_fold_change_proc)(data_gb, gene, order) for gene in tqdm(data.columns[:-1], desc="Calculating Log Fold Change", unit="gene")
    )
    
    result_df = pd.DataFrame(data=result, columns=['Ensembl_ID', 'logfc']).set_index('Ensembl_ID')
    return result_df

study = 'E-GEOD-111727'

adata = sc.datasets.ebi_expression_atlas(study)

# Basic filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Library size correction
sc.pp.normalize_total(adata, target_sum=1e4)

# Logarithmize the data
sc.pp.log1p(adata)

metadata = adata.obs
adata = adata.to_df()

label = metadata['Factor Value[disease]']

results = anova_genes_scipy(adata, label, cpu_count()-1)
results = qvalues(results)

logfc = log_fold_change(adata, label, 10)
results['logfc'] = logfc['logfc']

output_path = "../../data/single_cell"
disease = "HIV"

# Ensure the output directory exists
os.makedirs(output_path, exist_ok=True)

# Save the results DataFrame to a CSV file
results.to_csv(os.path.join(output_path, f'{disease}_{study}.csv'), index=True)