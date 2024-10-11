import pandas as pd
import numpy as np
import scanpy as sc
from src.util.pathway import parse_gmt

#study = 'E-ANND-3'
#data = sc.datasets.ebi_expression_atlas(study)

data = sc.read_mtx('data/E-ANND-3/E-ANND-3.aggregated_filtered_normalised_counts.mtx')
data = data.transpose()
data.var_names = pd.read_csv('data/E-ANND-3/E-ANND-3.aggregated_filtered_normalised_counts.mtx_rows', sep='\t', header=None).iloc[:,0]
data.obs = pd.read_csv('data/E-ANND-3/experimental_design.tsv', sep='\t', index_col=0)


groups = [*np.repeat(1,np.round(data.n_obs/2)),*np.repeat(2,data.n_obs - np.round(data.n_obs/2))]
data.obs['group'] = np.random.permutation(groups)

group1 = data.obs[data.obs['group'] == 1].index

kegg = parse_gmt('data/singlecell_kegg/KEGG_2022_entrez.gmt') #this parser is broken :(
pathway_genes = kegg.loc['hsa05218', 'genes']


delta = 10

data[group1,pathway_genes] = data[group1,pathway_genes] + delta
