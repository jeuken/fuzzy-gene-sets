import scanpy as sc
import os
import numpy as np
import pandas as pd
from scipy.stats import f_oneway
from joblib import Parallel, delayed, cpu_count
from qvalue import qvalues
from tqdm import tqdm  

# Define constants
STUDY = 'E-CURD-126'
KEGG_PATHWAY_NAME = 'hsa05170'
OUTPUT_DIR = "../../data/single_cell"
NEGATIVE_LABEL = "group2"
POSITIVE_LABEL = "group1"
DISEASE_NAME = "HIV"

def load_and_filter_data(study):
    adata = sc.datasets.ebi_expression_atlas(study)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"Number of cells after filtering: {adata.n_obs}")
    print(f"Number of genes after filtering: {adata.n_vars}")
    
    # Normalize and filter for normal samples
    sc.pp.normalize_total(adata, target_sum=1e4)
    adata = adata[adata.obs['Factor Value[disease]'] == 'normal']
    print("Disease status distribution after filtering:")
    print(adata.obs['Factor Value[disease]'].value_counts())
    
    return adata

def assign_groups(adata):
    np.random.seed(42)
    num_cells = adata.shape[0]
    adata.obs['group'] = np.random.choice([NEGATIVE_LABEL, POSITIVE_LABEL], size=num_cells)
    print("Group assignment distribution:")
    print(adata.obs['group'].value_counts())

def filter_pathway_genes(adata):
    kegg = pd.read_csv('/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/KEGG_2022_ensembl.csv')
    pathway_genes = kegg.loc[kegg['Pathway_Name'] == KEGG_PATHWAY_NAME, 'Ensembl_ID'].values
    genes_present = [gene for gene in pathway_genes if gene in adata.var_names]

    print(f"Number of genes in the pathway: {len(pathway_genes)}")
    print(f"Number of pathway genes present in the dataset: {len(genes_present)}")
    print("First 5 genes present in the dataset:")
    print(genes_present[:5])
    
    return genes_present

def modify_expression(adata, genes_present):
    group1_cells = adata.obs['group'] == POSITIVE_LABEL
    adata[group1_cells, genes_present].X *= 1.2  # Ensure this modification makes biological sense
    sc.pp.log1p(adata)
    print("Modified expression matrix for group1 (first 5 cells and first 5 genes):")
    print(adata[group1_cells, genes_present].X[:5, :5])

def perform_anova_and_logfc(adata, labels):
    adata.obs['label'] = labels
    metadata_filtered = adata.obs[adata.obs['label'].isin([NEGATIVE_LABEL, POSITIVE_LABEL])]
    adata_filtered = adata[metadata_filtered.index, :]
    label_filtered = metadata_filtered['label']
    
    n_jobs = cpu_count() - 1
    results = anova_genes_scipy(adata_filtered.to_df(), label_filtered, n_jobs)
    results = qvalues(results)
    
    logfc = log_fold_change(adata_filtered.to_df(), label_filtered, n_jobs, POSITIVE_LABEL)
    results['logfc'] = logfc['logfc']
    
    return results

def anova_proc_scipy(data_groupby, gene_name):
    groups = data_groupby[gene_name].apply(list)
    if any(len(set(group)) == 1 for group in groups):
        return (gene_name, np.nan, np.nan)
    
    f, p = f_oneway(*groups)
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
    mean1 = np.mean(unpacked[order[0]]) if order[0] == positive_label else np.mean(unpacked[order[1]])
    mean2 = np.mean(unpacked[order[1]]) if order[0] == positive_label else np.mean(unpacked[order[0]])
    fc = mean1 - mean2
    return (gene_name, fc)

def log_fold_change(data, label, n_jobs, positive_label):
    # Ensure the label is of categorical type
    data['label'] = label.astype('category')  # Convert to categorical
    data_gb = data.groupby('label', observed=False)
    order = data['label'].cat.categories  # Now this will work without raising an error

    result = Parallel(n_jobs=n_jobs)(
        delayed(log_fold_change_proc)(data_gb, gene, order, positive_label) for gene in tqdm(data.columns[:-1], desc="Calculating Log Fold Change", unit="gene")
    )

    result_df = pd.DataFrame(data=result, columns=['Ensembl_ID', 'logfc']).set_index('Ensembl_ID')
    return result_df

def save_results(results):
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    results.to_csv(os.path.join(OUTPUT_DIR, f'{DISEASE_NAME}_synthetic.csv'), index=True)
    print(f"Results saved to {os.path.join(OUTPUT_DIR, f'{DISEASE_NAME}_synthetic.csv')}")

def main():
    adata = load_and_filter_data(STUDY)
    assign_groups(adata)
    genes_present = filter_pathway_genes(adata)
    modify_expression(adata, genes_present)
    
    results = perform_anova_and_logfc(adata, adata.obs['group'])
    save_results(results)

if __name__ == "__main__":
    main()

