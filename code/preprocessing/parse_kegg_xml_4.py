# Import necessary libraries for parsing and data manipulation
import xmltodict  # To parse XML files into Python dictionaries.
import pandas as pd  # For data manipulation.
import networkx as nx  # For graph operations.
import matplotlib.pyplot as plt  # For plotting.

# Reading the KGML file and parsing it into a Python dictionary.
with open('/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/kgml/hsa05414.kgml', 'r') as file:
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
# Create an initial DataFrame for gene-to-gene relationships
gene_gene = pd.DataFrame(columns=['gene1', 'gene2', 'subtype'])

# Function to extract subtype
def extract_subtype(subtype):
    if isinstance(subtype, list) and len(subtype) > 0:
        return subtype[0]['@name'] if '@name' in subtype[0] else None
    elif isinstance(subtype, dict) and '@name' in subtype:
        return subtype['@name']
    return subtype   

# Handle all types of PPel relationships: Here please first filter for PPrel relations
for _, relation in relations.iterrows():
    entry1 = relation['@entry1']
    entry2 = relation['@entry2']
    
    # Check if entry1 is a gene and entry2 is a gene
    if entry1 in genes['@id'].values and entry2 in genes['@id'].values:
        gene1_name = id_to_name.get(entry1, entry1)
        gene2_name = id_to_name.get(entry2, entry2)
        
        # Check for subtype and extract name if available
        subtype = relation.get('subtype')  # Access subtype directly
        if isinstance(subtype, list) and len(subtype) > 0:
            subtype = subtype[0]['@name'] if '@name' in subtype[0] else None
        elif isinstance(subtype, dict) and '@name' in subtype:
            subtype = subtype['@name']
        
        gene_gene = pd.concat([gene_gene, pd.DataFrame([[gene1_name, gene2_name, subtype]],
                                                        columns=['gene1', 'gene2', 'subtype'])], ignore_index=True)

    elif entry1 in genes['@id'].values and entry2 in entries[entries['@type'] == 'group']['@id'].values:
        print(f"Found relation: {entry1} (gene) to {entry2} (group)")
        gene1_name = id_to_name.get(entry1, entry1)
        
        # Access the 'component' for the entry2 (which is a group)
        group_genes = entries.loc[entries['@id'] == entry2, 'component'].values[0]  # Access the first element
        print(group_genes)
        
        # Extract gene IDs from the group_genes
        group_gene_ids = [gene['@id'] for gene in group_genes if isinstance(gene, dict) and '@id' in gene]
        print(group_gene_ids)
        
        # Add relationships from gene to each gene in the group
        for group_gene_id in group_gene_ids:
            group_gene_name = id_to_name.get(group_gene_id, group_gene_id)
            gene_gene = pd.concat([gene_gene, pd.DataFrame([[gene1_name, group_gene_name, 'gene_to_group']],
                                                            columns=['gene1', 'gene2', 'subtype'])], ignore_index=True)

    # Check if entry1 is a group and entry2 is a gene
    elif entry1 in entries[entries['@type'] == 'group']['@id'].values and entry2 in genes['@id'].values:
        print(f"Found relation: {entry1} (group) to {entry2} (gene)")
        
        # Access the 'component' for the group (entry1)
        group_gene_ids = entries.loc[entries['@id'] == entry1, 'component'].values[0]  # Access the first element
        print(group_gene_ids)
        
        # Extract gene IDs from the group_gene_ids
        group_gene_ids = [gene['@id'] for gene in group_gene_ids if isinstance(gene, dict) and '@id' in gene]
        print(group_gene_ids)
        
        # Add relationships from each gene in the group to gene2
        for group_gene_id in group_gene_ids:
            group_gene_name = id_to_name.get(group_gene_id, group_gene_id)
            gene2_name = id_to_name.get(entry2, entry2)
            new_entries = pd.DataFrame([[group_gene_name, gene2_name, 'group_to_gene']],
                                       columns=['gene1', 'gene2', 'subtype'])
            gene_gene = pd.concat([gene_gene, new_entries], ignore_index=True)

    # Check if both entry1 and entry2 are groups
    elif entry1 in entries[entries['@type'] == 'group']['@id'].values and entry2 in entries[entries['@type'] == 'group']['@id'].values:
        print(f"Found relation: {entry1} (group) to {entry2} (group)")
        
        # Access the 'component' for the group (entry1)
        group1_gene_ids = entries.loc[entries['@id'] == entry1, 'component'].values[0]  # Access the first element
        group1_gene_ids = [gene['@id'] for gene in group1_gene_ids if isinstance(gene, dict) and '@id' in gene]
        print(f"Group1 genes: {group1_gene_ids}")
        
        # Access the 'component' for the group (entry2)
        group2_gene_ids = entries.loc[entries['@id'] == entry2, 'component'].values[0]  # Access the first element
        group2_gene_ids = [gene['@id'] for gene in group2_gene_ids if isinstance(gene, dict) and '@id' in gene]
        print(f"Group2 genes: {group2_gene_ids}")

        # Create relationships between all genes in group1 and group2
        for g1 in group1_gene_ids:
            for g2 in group2_gene_ids:
                g1_name = id_to_name.get(g1, g1)
                g2_name = id_to_name.get(g2, g2)
                new_entries = pd.DataFrame([[g1_name, g2_name, 'group_to_group']],
                                           columns=['gene1', 'gene2', 'subtype'])
                gene_gene = pd.concat([gene_gene, new_entries], ignore_index=True)
                print(f"Added group-to-group relation: {g1_name} -> {g2_name}")


# Handle all types of PCrel relationships: Here please filter for Pcrel relations
pcrel_rels = relations[relations['@type'] == 'PCrel']  # Select PCrel relations

for _, pcrel1 in pcrel_rels.iterrows():
    entry1 = pcrel1['@entry1']
    compound_id = pcrel1['@entry2']

    # Check if entry1 is a gene
    if entry1 in genes['@id'].values:
        gene1_name = id_to_name.get(entry1, entry1)
        
        # Check if the compound has associated genes or groups
        associated_entries = pcrel_rels[pcrel_rels['@entry1'] == compound_id]['@entry2'].tolist()
        
        for entry2 in associated_entries:
            if entry2 in genes['@id'].values:  # Entry2 is a gene
                gene2_name = id_to_name.get(entry2, entry2)
                indirect_rel = pd.DataFrame([[gene1_name, gene2_name, 'indirect_relation']],
                                             columns=['gene1', 'gene2', 'subtype'])
                gene_gene = pd.concat([gene_gene, indirect_rel], ignore_index=True)
                
            elif entry2 in entries[entries['@type'] == 'group']['@id'].values:  # Entry2 is a group
                group_genes = entries.loc[entries['@id'] == entry2, 'component'].values[0]
                group_gene_ids = [gene['@id'] for gene in group_genes if isinstance(gene, dict) and '@id' in gene]

                # Create relationships to each gene in the group
                for group_gene_id in group_gene_ids:
                    group_gene_name = id_to_name.get(group_gene_id, group_gene_id)
                    indirect_rel = pd.DataFrame([[gene1_name, group_gene_name, 'indirect_relation']],
                                                 columns=['gene1', 'gene2', 'subtype'])
                    gene_gene = pd.concat([gene_gene, indirect_rel], ignore_index=True)

    # Check if entry1 is a group
    elif entry1 in entries[entries['@type'] == 'group']['@id'].values:
        group_gene_ids = entries.loc[entries['@id'] == entry1, 'component'].values[0]
        group_gene_ids = [gene['@id'] for gene in group_gene_ids if isinstance(gene, dict) and '@id' in gene]

        # For each gene in the group, find indirect relations to the compound
        for group_gene_id in group_gene_ids:
            group_gene_name = id_to_name.get(group_gene_id, group_gene_id)
            # Check if the compound has associated genes or groups
            associated_entries = pcrel_rels[pcrel_rels['@entry1'] == compound_id]['@entry2'].tolist()
            
            for entry2 in associated_entries:
                if entry2 in genes['@id'].values:  # Entry2 is a gene
                    gene2_name = id_to_name.get(entry2, entry2)
                    indirect_rel = pd.DataFrame([[group_gene_name, gene2_name, 'indirect_relation']],
                                                 columns=['gene1', 'gene2', 'subtype'])
                    gene_gene = pd.concat([gene_gene, indirect_rel], ignore_index=True)
                    
                elif entry2 in entries[entries['@type'] == 'group']['@id'].values:  # Entry2 is a group
                    nested_group_genes = entries.loc[entries['@id'] == entry2, 'component'].values[0]
                    nested_group_gene_ids = [gene['@id'] for gene in nested_group_genes if isinstance(gene, dict) and '@id' in gene]

                    # Create relationships to each gene in the nested group
                    for nested_group_gene_id in nested_group_gene_ids:
                        nested_group_gene_name = id_to_name.get(nested_group_gene_id, nested_group_gene_id)
                        indirect_rel = pd.DataFrame([[group_gene_name, nested_group_gene_name, 'indirect_relation']],
                                                     columns=['gene1', 'gene2', 'subtype'])
                        gene_gene = pd.concat([gene_gene, indirect_rel], ignore_index=True)


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

# Split gene1 and gene2 values on spaces to ensure they are lists
gene_gene['gene1'] = gene_gene['gene1'].apply(lambda x: x.split() if isinstance(x, str) else [x])
gene_gene['gene2'] = gene_gene['gene2'].apply(lambda x: x.split() if isinstance(x, str) else x)

# Explode the gene_gene DataFrame to have individual gene entries
gene_gene_explode = gene_gene.explode('gene1').explode('gene2')

# Prepare the graph from the gene-to-gene relationships
G = nx.DiGraph()  # Create a directed graph

# Add edges based on gene-to-gene relationships
for _, row in gene_gene_explode.iterrows():
    G.add_edge(row['gene1'], row['gene2'])  # Add an edge from gene1 to gene2

# Set of all gene names
all_genes = set(genes['@name'])  # Set of all genes from the genes DataFrame
related_genes = set(gene_gene_explode['gene1']).union(set(gene_gene_explode['gene2']))  # Set of genes involved in relationships
isolated_genes = all_genes - related_genes  # Genes with no relationships

# Handle isolated genes
for gene in isolated_genes:
    G.add_node(gene)

# Calculate centrality metrics
centrality = nx.betweenness_centrality(G)  # Calculate betweenness centrality
degree_centrality = nx.degree_centrality(G)  # Calculate degree centrality

# Prepare to store results
centrality_df = pd.DataFrame.from_dict(centrality, orient='index', columns=['betweenness_centrality'])
degree_centrality_df = pd.DataFrame.from_dict(degree_centrality, orient='index', columns=['degree_centrality'])

# Save the gene-gene relationship DataFrames to separate CSV files
gene_gene.to_csv('../../data/pathways/KEGG/gene_gene_relationships.csv', index=False)
gene_gene_explode.to_csv('../../data/pathways/KEGG/gene_gene_relationships_exploded.csv', index=False)

# Plotting the directed graph
plt.figure(figsize=(12, 12))  # Set figure size
pos = nx.spring_layout(G)  # Set layout for nodes
nx.draw(G, pos, with_labels=True, node_size=2000, node_color='skyblue', font_size=10, font_color='black', font_weight='bold', arrows=True)

# Save the plot
plt.title('Gene Interaction Network')
plt.savefig('../../data/pathways/KEGG/gene_interaction_network.png', format='png')  # Save the plot as PNG
plt.show()  # Display the plot

