import os
import xmltodict  # To parse XML files into Python dictionaries.
import pandas as pd  # For data manipulation.
import networkx as nx  # For graph operations.
import matplotlib.pyplot as plt  # For plotting.
from joblib import Parallel, delayed  # For parallel processing
from tqdm import tqdm

# Function to process each KGML file
def process_kgml_file(kgml_file):
    # Attempt to parse the KGML file
    try:
        with open(kgml_file, 'r') as file:
            path_dict = xmltodict.parse(file.read())
    except Exception as e:
        print(f"Error processing file {kgml_file}: {e}")
        return pd.DataFrame()  # Return an empty DataFrame

    # Extract pathway name and description
    pathway_name = path_dict['pathway']['@name'].replace('path:', '')
    pathway_description = path_dict['pathway'].get('@title', '')  # Use get to handle missing description

    # Creating DataFrames for 'entry' nodes
    # Creating DataFrames for 'entry' nodes
    entries = pd.DataFrame.from_dict(path_dict['pathway']['entry'])  # Extract entries (genes)
    genes = entries[entries['@type'] == 'gene']  # Filter to include only gene entries
    genes.loc[:, '@name'] = genes['@name'].str.replace('hsa:', '', regex=False)  # Clean gene names

    # Initialize gene_gene DataFrame
    gene_gene = pd.DataFrame(columns=['gene1', 'gene2', 'subtype'])  # Initialize gene_gene DataFrame
    id_to_name = dict(zip(genes['@id'], genes['@name']))  # Map gene IDs to gene names
    relations = pd.DataFrame.from_dict(path_dict['pathway'].get('relation', []))  # Extract relations
    
    reactions = pd.DataFrame.from_dict(path_dict['pathway'].get('reaction', [])) 
    
    # Function to extract subtype
    def extract_subtype(subtype):
        if isinstance(subtype, list) and len(subtype) > 0:
            return subtype[0]['@name'] if '@name' in subtype[0] else None
        elif isinstance(subtype, dict) and '@name' in subtype:
            return subtype['@name']
        return subtype

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
        group_gene_list = []
        for _, group in groups.iterrows():
            gene_ids = [gene['@id'] for gene in group['component'] if isinstance(gene, dict) and '@id' in gene]
            for i in range(len(gene_ids)):
                for j in range(i + 1, len(gene_ids)):  # Avoid duplicate pairs
                    gene1_name = id_to_name.get(gene_ids[i], gene_ids[i])
                    gene2_name = id_to_name.get(gene_ids[j], gene_ids[j])
                    group_gene_list.append([gene1_name, gene2_name, 'group_association'])

        if group_gene_list:  # Create DataFrame if there are group relations
            gene_gene = pd.concat([gene_gene, pd.DataFrame(group_gene_list, columns=['gene1', 'gene2', 'subtype'])], ignore_index=True)

    gene_gene['gene1'] = gene_gene['gene1'].apply(lambda x: x.split() if isinstance(x, str) else [x])
    gene_gene['gene2'] = gene_gene['gene2'].apply(lambda x: x.split() if isinstance(x, str) else x)
    
    # Explode the gene_gene DataFrame to have individual gene entries
    gene_gene_explode = gene_gene.explode('gene1').explode('gene2').reset_index(drop=True)  # Explode both gene1 and gene2

    # Create directed graph
    G = nx.DiGraph()
    
    # Add edges from the exploded gene_gene DataFrame
    G.add_edges_from(gene_gene_explode[['gene1', 'gene2']].values)
    
    # Identify genes that have no relations
    all_genes = set(genes['@name'])  # All gene names from the KGML file
    connected_genes = set(G.nodes())  # Genes that are connected in the graph
    isolated_genes = all_genes - connected_genes  # Genes that are not connected
    
    # Add isolated genes as nodes in the graph without edges
    G.add_nodes_from(isolated_genes)
    
    # Calculate centrality metrics
    degree_dict = dict(G.degree())
    betweenness_dict = nx.betweenness_centrality(G)
    
    # Create DataFrame for results
    centrality_df = pd.DataFrame({
        'gene': list(degree_dict.keys()),
        'degree': list(degree_dict.values()),
        'betweenness': [betweenness_dict.get(gene, 0) for gene in degree_dict.keys()],  # Default to 0 for isolated nodes
        'pathway_name': pathway_name,  # Add pathway name
        'pathway_description': pathway_description  # Add pathway description
    })
    
    # Explode gene names if they contain multiple IDs
    centrality_df['gene'] = centrality_df['gene'].apply(lambda x: x.split())
    long_centrality_df = centrality_df.explode('gene').reset_index(drop=True)

    # Specify the directory for saving plots
    plot_directory = '../../data/pathways/KEGG/pathway_plots/'
    os.makedirs(plot_directory, exist_ok=True)
    
    # Derive file name for saving plots
    file_name = os.path.basename(kgml_file).replace('.kgml', '.png')
    
    # Plot the graph
    plt.figure(figsize=(15, 10))
    pos = nx.spring_layout(G, k=0.5)  # Set layout for the graph
    nx.draw_networkx_nodes(G, pos, node_size=500, node_color='lightblue')
    nx.draw_networkx_edges(G, pos, arrowstyle='-|>', arrowsize=10)
    
    # Modify labels to only include the first part of the gene name
    labels = {node: node.split()[0] for node in G.nodes()}  
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8, font_color='black')
    
    # Title the plot with the pathway name
    plt.title(f"Topology for {pathway_name}: {pathway_description}")
    plt.axis('off')  # Turn off the axis
    
    # Save the plot with the pathway name and ID
    plt.savefig(os.path.join(plot_directory, file_name), format='png')
    plt.close()  # Close the plot to free memory
    
    # Return long centrality DataFrame
    return long_centrality_df

# Use joblib for parallel processing
pathway_files = [os.path.join('../../data/pathways/KEGG/kgml/', file) for file in os.listdir('../../data/pathways/KEGG/kgml/') if file.endswith('.kgml')]
results = Parallel(n_jobs=-1)(delayed(process_kgml_file)(kgml_file) for kgml_file in tqdm(pathway_files))

# Combine all centrality DataFrames into one DataFrame
combined_centrality_df = pd.concat(results, ignore_index=True)
combined_centrality_df.to_csv('../../data/pathways/KEGG/kegg_centrality.csv', index=False)