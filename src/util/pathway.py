import pandas as pd

def parse_gmt(gmt_path, gene_id_dict_path = None):
    #gene_id_dict_path is the path to a file exported from biomart with only Gene Stable ID and Gene Name fields
    if gene_id_dict_path is not None:
        dict = pd.read_csv(gene_id_dict_path, sep='\t', index_col=1).dropna()
    gmt = pd.DataFrame(columns=['genes', 'url'])
    with open(gmt_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.split('\t')
            gmt.loc[line[0],'url'] = line[1]
            if gene_id_dict_path is not None:
                genes = dict[dict.index.isin(line[2:])].values.T.tolist()[0]
            else:
                genes = line[2:]
            gmt.loc[line[0],'genes'] = genes
    return gmt
