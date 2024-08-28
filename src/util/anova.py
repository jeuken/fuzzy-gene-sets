from scipy.stats import f_oneway
from joblib import Parallel, delayed, cpu_count

def anova_proc_scipy(data_groupby, gene_name):
    f, p = f_oneway(*data_groupby[gene_name].apply(list))
    return (gene_name, p)

def anova_genes_scipy(data, label, n_jobs):
    data['label'] = label 
    data_gb = data.groupby('label')
    result = Parallel(n_jobs=n_jobs)(delayed(anova_proc_scipy)(data_gb, gene) for gene in data.columns[:-1])
    result_df = pd.DataFrame(data = result, columns=['gene', 'p']).set_index('gene')
    return result_df
