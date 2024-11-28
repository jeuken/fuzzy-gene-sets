import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import f_oneway
from joblib import Parallel, delayed, cpu_count
from src.util.qvalue import qvalues



def anova_proc_scipy(data_groupby, gene_name):
    f, p = f_oneway(*data_groupby[gene_name].apply(list))
    return (gene_name, p, f)

def anova_genes_scipy(data, label, n_jobs):
    data['label'] = label 
    data_gb = data.groupby('label')
    result = Parallel(n_jobs=n_jobs)(delayed(anova_proc_scipy)(data_gb, gene) for gene in data.columns[:-1])
    result_df = pd.DataFrame(data = result, columns=['gene', 'p', 'f']).set_index('gene')
    return result_df

def log_fold_change_proc(data_groupby, gene_name, order):
    unpacked = data_groupby[gene_name].apply(list)
    mean1 = np.mean(unpacked[order[0]])
    mean2 = np.mean(unpacked[order[1]])
    fc = mean1 - mean2
    return (gene_name, fc)

def log_fold_change(data, label, n_jobs):
    data['label'] = label 
    data_gb = data.groupby('label')
    order = data_gb[data.columns[0]].apply(list).index
    result = Parallel(n_jobs=n_jobs)(delayed(log_fold_change_proc)(data_gb, gene, order) for gene in data.columns[:-1])
    result_df = pd.DataFrame(data = result, columns=['gene', 'logfc',]).set_index('gene')
    return result_df


def diff_proc(adata, i, label):
    exp_df = pd.DataFrame(adata.X.getcol(i).todense())
    exp_df['label'] = label.values
    exp_groupby = exp_df.groupby('label', observed=True)
    f, p = f_oneway(*exp_groupby[0].apply(list))
    if i%100 == 0:
        print(i)
    return (i, p, f)

def diff_par(adata, label, n_jobs):
    n_genes = adata.X.shape[1]
    result = Parallel(n_jobs=n_jobs)(delayed(diff_proc)(adata, i, label) for i in range(n_genes))
    result_df = pd.DataFrame(result, columns=['gene_i', 'p', 'f']).set_index('gene_i')
    return result_df


studies = pd.read_csv('./data/singlecell_kegg/sc_kegg_target_disease.tsv', sep = '\t', index_col = 0)

# study = 'E-GEOD-76312'
# study = studies.index[15]

for study in studies.index:
    print(study)

    adata = sc.datasets.ebi_expression_atlas(study)

    # Basic filtering
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    # Library size correction
    sc.pp.normalize_total(adata, target_sum=1e4)

    # Logarithmize the data
    sc.pp.log1p(adata)

    metadata = adata.obs
    # adata = adata.to_df()

    label = metadata['Factor Value[disease]']

    # results = anova_genes_scipy(adata, label, cpu_count()-1)
    results = diff_par(adata, label, cpu_count()-1)
    results = qvalues(results)
    results['gene_name'] = adata.var.iloc[results.index].index
    results = results.set_index('gene_name')

    results.to_csv(study + '_result.csv')

