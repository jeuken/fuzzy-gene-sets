# Import necessary libraries for parsing and data manipulation
import xmltodict  # To parse XML files into Python dictionaries.
import pandas as pd  # For data manipulation.
import networkx as nx  # For graph operations.
import matplotlib.pyplot as plt  # For plotting.

# Reading the KGML file and parsing it into a Python dictionary.
with open('/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/kgml/hsa05414.kgml', 'r') as file:
    path_dict = xmltodict.parse(file.read())

# Creating DataFrames for 'entry', 'reaction', and 'relation' nodes.
entries = pd.DataFrame.from_dict(path_dict['pathway']['entry'])  # Extract entries (genes).

# Initialize DataFrames for reactions and relations
reactions = pd.DataFrame.from_dict(path_dict['pathway'].get('reaction', []))  # Extract reactions, default to empty list
relations = pd.DataFrame.from_dict(path_dict['pathway'].get('relation', []))  # Extract relations, default to empty list

# Filtering to include only gene entries
genes = entries[entries['@type'] == 'gene']
compounds = entries[entries['@type'] == 'compound']
enzymes = entries[entries['@type'] == 'enzyme']
groups = entries[entries['@type'] == 'group']

# Clean gene names by removing 'hsa:' prefix
genes['@name'] = genes['@name'].str.replace('hsa:', '', regex=False)

# DataFrame to store gene-related information: names, substrates, and products
gene_product_substrate = pd.DataFrame(columns=['name', 'substrate', 'product'])

# Extracting reaction relationships for each gene
for gene_index in genes.index:
    gene_name = genes.loc[gene_index, '@name']
    gene_product_substrate.loc[gene_index, 'name'] = gene_name

    # Getting the associated reaction for this gene, if it exists
    reaction = genes.loc[gene_index, '@reaction'] if '@reaction' in genes.columns else None

    # Extract substrates (input molecules) for the reaction
    if reaction and not reactions.empty:
        substrates = reactions.loc[reactions['@name'] == reaction, 'substrate'].iloc[0]
        substrate_ids = [s['@id'] for s in substrates] if isinstance(substrates, list) else [substrates['@id']]
        gene_product_substrate.loc[gene_index, 'substrate'] = substrate_ids
    else:
        gene_product_substrate.loc[gene_index, 'substrate'] = []

    # Extract products (output molecules) for the reaction
    if reaction and not reactions.empty:
        products = reactions.loc[reactions['@name'] == reaction, 'product'].iloc[0]
        product_ids = [p['@id'] for p in products] if isinstance(products, list) else [products['@id']]
        gene_product_substrate.loc[gene_index, 'product'] = product_ids
    else:
        gene_product_substrate.loc[gene_index, 'product'] = []

# DataFrame to store gene-to-gene relationships based on substrates/products
gene_gene_list = []  # Use a list to collect gene relationships

# Establish direct connections based on substrate-product relationships
for gene_index1 in gene_product_substrate.index:
    products = gene_product_substrate.loc[gene_index1, 'product']
    
    for product in products:
        for gene_index2 in gene_product_substrate.index:
            if product in gene_product_substrate.loc[gene_index2, 'substrate']:
                gene_gene_list.append([
                    gene_product_substrate.loc[gene_index1, 'name'],
                    gene_product_substrate.loc[gene_index2, 'name'],
                    'substrate-product association'
                ])

# Convert list to DataFrame
gene_gene = pd.DataFrame(gene_gene_list, columns=['gene1', 'gene2', 'subtype'])
print("gene_gene_before:",gene_gene)
# Map gene IDs to gene names
id_to_name = dict(zip(genes['@id'], genes['@name']))

# Filter relations to include only gene pairs
if not relations.empty:
    # Get gene IDs
    gene_ids = set(genes['@id'])
    
    # Filter relations to include only gene pairs
    filtered_relations = relations[
        relations['@entry1'].isin(gene_ids) & relations['@entry2'].isin(gene_ids)
    ].copy()
    
    # Replace gene IDs with gene names in relations
    filtered_relations['gene1'] = filtered_relations['@entry1'].map(id_to_name)
    filtered_relations['gene2'] = filtered_relations['@entry2'].map(id_to_name)
    
    # Extract the '@name' from the 'subtype', using the first one if it's a list of dictionaries
    def extract_subtype(subtype):
        if isinstance(subtype, list) and len(subtype) > 0:
            return subtype[0]['@name'] if '@name' in subtype[0] else None
        elif isinstance(subtype, dict) and '@name' in subtype:
            return subtype['@name']
        else:
            return subtype
    
    filtered_relations['subtype'] = filtered_relations['subtype'].apply(extract_subtype)
    
    # Add filtered relations to gene_gene DataFrame
    gene_gene = pd.concat([gene_gene, filtered_relations[['gene1', 'gene2', 'subtype']]], ignore_index=True)

    # Handle PCrel to find indirect relations
    pcrel_rels = relations[relations['@type'] == 'PCrel']  # Select PCrel relations

    # Map compounds to associated genes
    indirect_rels = []  # Use a list to collect indirect relations
    # Handle PCrel to find indirect relations
    pcrel_rels = relations[relations['@type'] == 'PCrel']  # Select PCrel relations
    
    # Map compounds to associated genes
    indirect_rels = []  # Use a list to collect indirect relations
    for _, pcrel1 in pcrel_rels.iterrows():
        if pcrel1['@entry1'] in gene_ids:
            gene1 = pcrel1['@entry1']
            compound_id = pcrel1['@entry2']
            # Find associated genes via compounds
            associated_genes = pcrel_rels[pcrel_rels['@entry1'] == compound_id]['@entry2'].tolist()
            for gene2 in associated_genes:
                # Map gene IDs to gene names
                gene1_name = id_to_name.get(gene1, gene1)  # Get gene name from ID
                gene2_name = id_to_name.get(gene2, gene2)  # Get gene name from ID
                indirect_rels.append([gene1_name, gene2_name, 'indirect_relation'])
    
    # Add indirect relationships to the DataFrame
    if indirect_rels:  # Only create DataFrame if there are indirect relations
        gene_gene = pd.concat([gene_gene, pd.DataFrame(indirect_rels, columns=['gene1', 'gene2', 'subtype'])], ignore_index=True)

# New part: Add group gene associations
if not groups.empty:
    group_gene_list = []  # Use a list to collect group gene associations
    for _, group in groups.iterrows():
        group_genes = group['component']

        # Extract IDs from the component dictionaries
        gene_ids = [gene['@id'] for gene in group_genes if isinstance(gene, dict) and '@id' in gene]

        # Create pairs from gene IDs and add to group_gene_list
        for i in range(len(gene_ids)):
            for j in range(i + 1, len(gene_ids)):  # To avoid duplicate pairs
                gene1_name = id_to_name.get(gene_ids[i], gene_ids[i])  # Get gene name from ID
                gene2_name = id_to_name.get(gene_ids[j], gene_ids[j])  # Get gene name from ID
                group_gene_list.append([gene1_name, gene2_name, 'group_association'])  # Add pair with subtype

    if group_gene_list:  # Only create DataFrame if there are group relations
        gene_gene = pd.concat([gene_gene, pd.DataFrame(group_gene_list, columns=['gene1', 'gene2', 'subtype'])], ignore_index=True)
    
    # Add bidirectional relations for specified subtypes
    bidirectional_subtypes = ['binding/association', 'group_association']
    for _, row in gene_gene.iterrows():
        gene1 = row['gene1']
        gene2 = row['gene2']
        subtype = row['subtype']
        
        # If subtype is bidirectional, add the reverse relation
        if subtype in bidirectional_subtypes:
            gene_gene = pd.concat([gene_gene, pd.DataFrame([[gene2, gene1, subtype]], columns=['gene1', 'gene2', 'subtype'])], ignore_index=True)

# Initialize directed graph
G = nx.DiGraph()

# Add edges from gene_gene DataFrame
G.add_edges_from([(row['gene1'], row['gene2']) for _, row in gene_gene.iterrows()])

# Identify genes that have no relations or reactions
all_genes = set(genes['@name'])  # All gene names
connected_genes = set(G.nodes())  # Genes that are connected in the graph
isolated_genes = all_genes - connected_genes  # Genes that are not connected

# Add isolated genes as nodes in the graph without edges
G.add_nodes_from(isolated_genes)

# Calculate node degree and betweenness centrality
degree_dict = dict(G.degree())  # Node degree
betweenness_dict = nx.betweenness_centrality(G)  # Betweenness centrality

# Create DataFrame for results
centrality_df = pd.DataFrame({
    'gene': list(degree_dict.keys()),
    'degree': list(degree_dict.values()),
    'betweenness': [betweenness_dict[gene] for gene in degree_dict.keys()]
})


# Extract the pathway name from path_dict
pathway_name = path_dict['pathway']['@title']  # Assuming '@title' contains the pathway name
file_name = f"{pathway_name.replace(' ', '_')}_graph.png"

# Plot the graph
plt.figure(figsize=(15, 10))
pos = nx.spring_layout(G, k=0.5)  # Set the layout for the graph
nx.draw_networkx_nodes(G, pos, node_size=500, node_color='lightblue')
nx.draw_networkx_edges(G, pos, arrowstyle='-|>', arrowsize=10)

# Modify labels to only include the first number
labels = {node: node.split()[0] for node in G.nodes()}  # Keep only the first part of the gene name
nx.draw_networkx_labels(G, pos, labels=labels, font_size=8, font_color='black')

# Title the plot with the pathway name
plt.title(f"Gene Interaction Network from {pathway_name}")
plt.axis('off')  # Turn off the axis

# Save the plot with the pathway name and improved resolution
plt.savefig(f"../../data/pathways/KEGG/pathway_plots/{file_name}", dpi=300)  # Increased DPI for higher resolution

# Show plot
plt.show()




