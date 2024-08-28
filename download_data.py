import scanpy as sc
import pandas as pd

datasets = pd.read_csv('data/singlecell_kegg/sc_kegg_target_disease.tsv', sep = '\t', index_col=0)

for study in datasets.index:
    adata = sc.datasets.ebi_expression_atlas(study)


