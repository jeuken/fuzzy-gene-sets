import xmltodict  # Importing the xmltodict library to parse XML files into Python dictionaries.
import pandas as pd  # Importing pandas for data manipulation.
import networkx as nx  # Importing networkx for graph operations.
import matplotlib.pyplot as plt  # Importing matplotlib for plotting.

# Reading the KGML (KEGG Markup Language) file and parsing its XML content into a Python dictionary.
with open('/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/kgml/hsa04024.kgml', 'r') as file:
    path_dict = xmltodict.parse(file.read())

# Creating dataframes for 'entry', 'reaction', and 'relation' nodes from the pathway dictionary.
entries = pd.DataFrame.from_dict(path_dict['pathway']['entry'])  # Extracting the 'entry' elements of the pathway (contains genes).

# Check if reactions and relations exist
reactions = pd.DataFrame()  # Initialize an empty DataFrame for reactions
relations = pd.DataFrame()  # Initialize an empty DataFrame for relations

if 'reaction' in path_dict['pathway']:
    reactions = pd.DataFrame.from_dict(path_dict['pathway']['reaction'])  # Extracting the 'reaction' elements of the pathway.
if 'relation' in path_dict['pathway']:
    relations = pd.DataFrame.from_dict(path_dict['pathway']['relation'])  # Extracting the 'relation' elements of the pathway.

# Filtering entries to only include those that are genes (entries where '@type' is 'gene').
genes = entries[entries['@type'] == 'gene']

# Initializing an empty dataframe to store gene-related information: gene names, substrates, and products.
gene_product_substrate = pd.DataFrame(columns=['name', 'substrate', 'product'])

# Iterating through the genes and extracting their reaction relationships.
for gene_index in genes.index:
    # Storing the gene name in the 'name' column.
    gene_name = genes.loc[gene_index, '@name']
    gene_product_substrate.loc[gene_index, 'name'] = gene_name
    
    # Getting the reaction associated with this gene, if it exists.
    reaction = genes.loc[gene_index, '@reaction'] if '@reaction' in genes.columns else None
    
    # Extracting substrates (input molecules) for the reaction, if reactions exist.
    if reaction and not reactions.empty:
        substrates = reactions[reactions['@name'] == reaction]['substrate']
        if not substrates.empty:
            substrates = substrates.iloc[0]
            if isinstance(substrates, list):
                substrate_ids = [s['@id'] for s in substrates]
            else:
                substrate_ids = [substrates['@id']]
            gene_product_substrate.loc[gene_index, 'substrate'] = substrate_ids
        else:
            gene_product_substrate.loc[gene_index, 'substrate'] = []
    else:
        gene_product_substrate.loc[gene_index, 'substrate'] = []

    # Extracting products (output molecules) for the reaction, if reactions exist.
    if reaction and not reactions.empty:
        products = reactions[reactions['@name'] == reaction]['product']
        if not products.empty:
            products = products.iloc[0]
            if isinstance(products, list):
                product_ids = [p['@id'] for p in products]
            else:
                product_ids = [products['@id']]
            gene_product_substrate.loc[gene_index, 'product'] = product_ids
        else:
            gene_product_substrate.loc[gene_index, 'product'] = []
    else:
        gene_product_substrate.loc[gene_index, 'product'] = []

# Initializing an empty dataframe to store gene-to-gene relationships based on substrates and products.
gene_gene = pd.DataFrame(columns=['gene1', 'gene2'])

idx = 0  # Initialize the index for gene-gene relationships

# Iterating through the genes to establish connections between them based on substrates and products.
for gene_index1 in gene_product_substrate.index:
    # Retrieving the products of the first gene.
    products = gene_product_substrate.loc[gene_index1, 'product']
    
    # For each product of the first gene, check if it is a substrate for another gene.
    for product in products:
        for gene_index2 in gene_product_substrate.index:
            if product in gene_product_substrate.loc[gene_index2, 'substrate']:
                # If a match is found, store the gene pair (gene1 -> gene2) in the 'gene_gene' dataframe.
                link = [gene_product_substrate.loc[gene_index1, 'name'], gene_product_substrate.loc[gene_index2, 'name']]
                gene_gene.loc[idx] = link  # Add the gene pair to the dataframe.
                idx += 1  # Increment the index for the next gene-gene relationship.
                
# Create a dictionary to map gene IDs to gene names
id_to_name = dict(zip(genes['@id'], genes['@name']))

# Filter relations to only include pairs where both @entry1 and @entry2 are in the gene_ids set, if relations exist.
if not relations.empty:
    gene_ids = set(genes['@id'])  # Get the set of gene IDs
    filtered_relations = relations[
        relations['@entry1'].isin(gene_ids) & relations['@entry2'].isin(gene_ids)
    ].copy() 

    # Replace gene IDs with gene names in the relations 
    filtered_relations.loc[:, 'gene1'] = filtered_relations['@entry1'].map(id_to_name)
    filtered_relations.loc[:, 'gene2'] = filtered_relations['@entry2'].map(id_to_name)

    # Use pd.concat() to combine with the existing gene_gene DataFrame
    gene_gene = pd.concat(
        [gene_gene, filtered_relations[['gene1', 'gene2']]],
        ignore_index=True
    )


# Initialize a directed graph
G = nx.DiGraph()

# Add edges from gene_gene DataFrame
for _, row in gene_gene.iterrows():
    G.add_edge(row['gene1'], row['gene2'])

# Calculate node degree and shortest path betweenness for each gene
degree_dict = dict(G.degree())  # Node degree
betweenness_dict = nx.betweenness_centrality(G)  # Betweenness centrality

# Create a DataFrame to store the results
centrality_df = pd.DataFrame({
    'gene': list(degree_dict.keys()),
    'degree': list(degree_dict.values()),
    'betweenness': [betweenness_dict[gene] for gene in degree_dict.keys()]
})

# Display the table of gene degree and betweenness
print(centrality_df)

# Plot the graph
plt.figure(figsize=(15, 15))
pos = nx.spring_layout(G, seed=42)  # Position the nodes using the spring layout
nx.draw(G, pos, with_labels=True, node_size=3000, node_color="lightblue", font_size=10, font_weight="bold", edge_color="gray", arrows=True)
plt.title("Gene-Gene Interaction Network", size=15)
plt.show()

# Extract '@name' from 'graphics' and split into a list of names
genes['names'] = genes['graphics'].apply(
    lambda x: x.get('@name', '').split(', ') if isinstance(x, dict) and '@name' in x else []
)

# Flatten the list of names into separate rows if needed
names_expanded = genes[['names']].explode('names').reset_index(drop=True)

print(genes[['names']])
print(names_expanded)