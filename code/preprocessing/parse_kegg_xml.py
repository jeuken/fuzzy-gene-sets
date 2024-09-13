import xmltodict
import pandas as pd

with open('../../data/hsa00010.kgml', 'r') as file:
    path_dict = xmltodict.parse(file.read())


entries = pd.DataFrame.from_dict(path_dict['pathway']['entry'])
reactions = pd.DataFrame.from_dict(path_dict['pathway']['reaction'])

genes = entries[entries['@type'] == 'gene']

gene_product_substrate = pd.DataFrame(columns=['name','substrate','product'])
for gene in genes.index:
    gene_product_substrate.loc[gene,'name'] = genes.loc[gene,'@name']
    reaction = genes.loc[gene,'@reaction']
    substrates = reactions[reactions['@name'] == reaction]['substrate']
    if isinstance(substrates.iloc[0], list):
        ids = []
        for substrate in substrates.iloc[0]:
            ids.append(substrate['@id'])
    else:
        ids = [substrates.iloc[0]['@id']]
    gene_product_substrate.loc[gene,'substrate'] = ids
    products = reactions[reactions['@name'] == reaction]['product']
    if isinstance(products.iloc[0], list):
        ids = []
        for product in products.iloc[0]:
            ids.append(product['@id'])
    else:
        ids = [products.iloc[0]['@id']]
    gene_product_substrate.loc[gene,'product'] = ids


gene_gene = pd.DataFrame(columns=['gene1','gene2'])
idx = 0
for gene_index1 in gene_product_substrate.index:
    products = gene_product_substrate.loc[gene_index1,'product']
    for product in products:
        for gene_index2 in gene_product_substrate.index:
            if product in gene_product_substrate.loc[gene_index2,'substrate']:
                link = [gene_product_substrate.loc[gene_index1, 'name'], gene_product_substrate.loc[gene_index2, 'name']]
                gene_gene.loc[idx] = link
                idx +=1

