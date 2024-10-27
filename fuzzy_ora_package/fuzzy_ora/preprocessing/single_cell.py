import scanpy as sc
import os
import pandas as pd
import numpy as np
from scipy.stats import f_oneway
from joblib import Parallel, delayed, cpu_count
from qvalue import qvalues
from tqdm import tqdm
import argparse


def anova_proc_scipy(data_groupby, gene_name):
    """Perform ANOVA on a given gene."""
    f, p = f_oneway(*data_groupby[gene_name].apply(list))
    return gene_name, p, f

def anova_genes_scipy(data, labels, n_jobs):
    """Run ANOVA for all genes in the dataset."""
    data['label'] = labels
    data_grouped = data.groupby('label', observed=False)
    
    results = Parallel(n_jobs=n_jobs)(
        delayed(anova_proc_scipy)(data_grouped, gene) for gene in tqdm(data.columns[:-1], desc="Running ANOVA", unit="gene")
    )
    
    return pd.DataFrame(results, columns=['Ensembl_ID', 'p', 'f']).set_index('Ensembl_ID')

def log_fold_change_proc(data_groupby, gene_name, order, positive_label):
    """Calculate log fold change for a given gene."""
    unpacked = data_groupby[gene_name].apply(list)
    mean1, mean2 = (np.mean(unpacked[order[0]]), np.mean(unpacked[order[1]])) if order[0] == positive_label else (np.mean(unpacked[order[1]]), np.mean(unpacked[order[0]]))
    return gene_name, mean1 - mean2

def log_fold_change(data, labels, n_jobs, positive_label):
    """Compute log fold changes for all genes."""
    data['label'] = labels
    data_grouped = data.groupby('label', observed=False)
    order = labels.cat.categories

    results = Parallel(n_jobs=n_jobs)(
        delayed(log_fold_change_proc)(data_grouped, gene, order, positive_label) for gene in tqdm(data.columns[:-1], desc="Calculating Log Fold Change", unit="gene")
    )
    
    return pd.DataFrame(results, columns=['Ensembl_ID', 'logfc']).set_index('Ensembl_ID')

def sc_exp_data_preprocessing(study, negative_label, positive_label, condition, output_path):
    """Main function to process expression data, run ANOVA, and calculate log fold change."""
    # Load dataset from Expression Atlas
    adata = sc.datasets.ebi_expression_atlas(study)

    # Basic filtering
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Library size correction and log transformation
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Get metadata and convert to DataFrame
    metadata = adata.obs
    adata_df = adata.to_df()

    # Ensure that label is categorical
    metadata['Factor Value[disease]'] = metadata['Factor Value[disease]'].astype('category')

    # Filter metadata and data based on selected labels
    metadata_filtered = metadata[metadata['Factor Value[disease]'].isin([negative_label, positive_label])]
    adata_filtered = adata_df.loc[metadata_filtered.index]
    label_filtered = metadata_filtered['Factor Value[disease]']

    # Ensure data contains genes before proceeding
    if adata_filtered.empty:
        raise ValueError("Filtered dataset contains no genes.")

    # Run ANOVA and log fold change calculations
    n_jobs = cpu_count() - 1
    anova_results = anova_genes_scipy(adata_filtered, label_filtered, n_jobs)
    qvalue_corrected = qvalues(anova_results['p'])

    logfc_results = log_fold_change(adata_filtered, label_filtered, n_jobs, positive_label)
    anova_results['logfc'] = logfc_results['logfc']
    anova_results['qvalue'] = qvalue_corrected

    # Ensure the output directory exists
    os.makedirs(output_path, exist_ok=True)

    # Save results to CSV
    output_file = os.path.join(output_path, f'{condition}_{study}.csv')
    anova_results.to_csv(output_file, index=True)

    print(f"Results saved to {output_file}")

def sc_exp_data_preprocessing_main():
    parser = argparse.ArgumentParser(description="Process single-cell data and perform ANOVA and log fold change calculations.")
    parser.add_argument('--study', type=str, required=True, help="Study ID to load from Expression Atlas (e.g., E-GEOD-111727)")
    parser.add_argument('--negative_label', type=str, required=True, help="Label for negative group (e.g., control)")
    parser.add_argument('--positive_label', type=str, required=True, help="Label for positive group (e.g., disease)")
    parser.add_argument('--condition', type=str, required=True, help="Short disease name for output files (e.g., HIV)")
    parser.add_argument('--output_path', type=str, default="../../data/single_cell", help="Path to save results")

    args = parser.parse_args()

    # Run main function with provided arguments
    sc_exp_data_preprocessing(args.study, args.negative_label, args.positive_label, args.condition, args.output_path)

if __name__ == "__main__":
    #sc_exp_data_preprocessing_main()
    study = 'E-GEOD-111727'
    condition = 'HIV'
    negative_label = "normal"
    positive_label = "HIV"
    output_path = "../../data/single_cell"
    
    
    
    