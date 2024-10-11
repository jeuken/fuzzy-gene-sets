## Import necessary libraries for parsing and data manipulation
import xmltodict  # To parse XML files into Python dictionaries.
import pandas as pd  # For data manipulation.
import networkx as nx  # For graph operations.
import matplotlib.pyplot as plt  # For plotting.

# Reading the KGML file and parsing it into a Python dictionary.
with open('/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/kgml/hsa00020.kgml', 'r') as file:
    path_dict = xmltodict.parse(file.read())

# Creating DataFrames for 'entry' nodes.
entries = pd.DataFrame.from_dict(path_dict['pathway']['entry'])  # Extract entries (genes).

# Filtering to include only gene entries
genes = entries[entries['@type'] == 'gene']

# Clean gene names by removing 'hsa:' prefix
genes.loc[:, '@name'] = genes['@name'].str.replace('hsa:', '', regex=False)

# DataFrame to store gene-to-gene relationships based on relations
gene_gene_list = []  # Use a list to collect gene relationships

# Map gene IDs to gene names
id_to_name = dict(zip(genes['@id'], genes['@name']))

# Initialize DataFrames for relations
relations = pd.DataFrame.from_dict(path_dict['pathway'].get('relation', []))  # Extract relations, default to empty list
reactions = pd.DataFrame.from_dict(path_dict['pathway'].get('reaction', [])) 
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
    gene_gene = filtered_relations[['gene1', 'gene2', 'subtype']]

    # Handle PCrel to find indirect relations
    pcrel_rels = relations[relations['@type'] == 'PCrel']  # Select PCrel relations
    for _, pcrel1 in pcrel_rels.iterrows():
        if pcrel1['@entry1'] in gene_ids:
            gene1 = pcrel1['@entry1']
            compound_id = pcrel1['@entry2']
            associated_genes = pcrel_rels[pcrel_rels['@entry1'] == compound_id]['@entry2'].tolist()
            indirect_rels = [
                [id_to_name.get(gene1, gene1), id_to_name.get(gene2, gene2), 'indirect_relation']
                for gene2 in associated_genes
            ]
            gene_gene = pd.concat([gene_gene, pd.DataFrame(indirect_rels, columns=['gene1', 'gene2', 'subtype'])], ignore_index=True)

# New part: Add group gene associations
groups = entries[entries['@type'] == 'group']  # Extract groups
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


# Assuming gene_gene is already created with gene1 and gene2 columns
# Assuming gene_gene is already created with gene1 and gene2 columns
# Split gene1 values on spaces to ensure they are lists
gene_gene['gene1'] = gene_gene['gene1'].apply(lambda x: x.split() if isinstance(x, str) else [x])

# Split gene2 values on spaces and keep them as lists
gene_gene['gene2'] = gene_gene['gene2'].apply(lambda x: x.split() if isinstance(x, str) else x)

# Explode the gene_gene DataFrame to have individual gene entries
gene_gene_explode = gene_gene.explode('gene1').explode('gene2').reset_index(drop=True)  # Explode both gene1 and gene2

# Now gene_gene_explode will have each pair of genes as separate rows


# Initialize directed graph
G = nx.DiGraph()

# Add edges from gene_gene DataFrame
G.add_edges_from([(row['gene1'], row['gene2']) for _, row in gene_gene_explode.iterrows()])

# Identify genes that have no relations
all_genes = set(genes['@name'])  # All gene names
connected_genes = set(G.nodes())  # Genes that are connected in the graph
isolated_genes = all_genes - connected_genes  # Genes that are not connected

# Add isolated genes as nodes in the graph without edges
G.add_nodes_from(isolated_genes)


# Extract the pathway name from path_dict
pathway_name = path_dict['pathway']['@title']  # Assuming '@title' contains the pathway name
pathway_id = path_dict['pathway']['@name'].replace('path:', '')
file_name = f"{pathway_name.replace(' ', '_')}_graph.png"

# Calculate node degree and betweenness centrality
degree_dict = dict(G.degree())  # Node degree
betweenness_dict = nx.betweenness_centrality(G)  # Betweenness centrality

# Create DataFrame for results
centrality_df = pd.DataFrame({
    'gene': list(degree_dict.keys()),
    'degree': list(degree_dict.values()),
    'betweenness': [betweenness_dict[gene] for gene in degree_dict.keys()]
})

# If any gene names have multiple IDs, split them and explode to separate rows
centrality_df['gene'] = centrality_df['gene'].apply(lambda x: x.split())  # Split genes into lists
long_centrality_df = centrality_df.explode('gene').reset_index(drop=True)  # Explode into separate rows
# Add Pathway_Name and Description columns to long_centrality_df
long_centrality_df['Pathway_Name'] = pathway_id  # Add Pathway_Name column
long_centrality_df['Description'] = pathway_name  # Add Description column

# Display the long centrality DataFrame
print(long_centrality_df)



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
