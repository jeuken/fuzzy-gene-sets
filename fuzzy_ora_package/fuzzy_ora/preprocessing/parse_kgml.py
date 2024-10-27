import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse
import xmltodict  # For parsing XML data into Python dictionaries
import networkx as nx  # For creating and manipulating complex networks
from joblib import Parallel, delayed  # For parallel computing
from tqdm import tqdm  # For displaying progress bars

def extract_relation_subtype(subtype):
    """Extracts the relation subtype from the subtype column"""
    if isinstance(subtype, list) and subtype:
        return subtype[0].get('@name')  # Get name from the first item in a non-empty list
    elif isinstance(subtype, dict):
        return subtype.get('@name')  # Get name from a dictionary
    return subtype  # Return subtype if neither condition is met

def extract_group_genes(group, gene_id_name_map):
    """Extracts gene names from a group, mapping IDs to names using the provided dictionary."""
    components = group['component'] if 'component' in group else []  # Get components if present
    gene_ids = [comp['@id'] for comp in components if isinstance(comp, dict) and '@id' in comp]  # Extract gene IDs
    return ' '.join([gene_id_name_map.get(gene_id, gene_id) for gene_id in gene_ids])  # Map IDs to names and join


def get_nodes_and_edges(path_dict):
    """
    Extract nodes and edges from the provided path dictionary.

    Parameters:
    - path_dict (dict): Dictionary containing pathway data.

    Returns:
    - nodes (DataFrame): DataFrame of nodes with their attributes.
    - edges (DataFrame): DataFrame of edges with source and target relationships.
    """
    # Extract pathway name and description
    pathway_name = path_dict['pathway']['@name'].replace('path:', '')  # Clean the pathway name
    pathway_description = path_dict['pathway'].get('@title', '')  # Get the pathway title, default to empty string

    # Create a DataFrame from pathway entries
    entries = pd.DataFrame.from_dict(path_dict['pathway']['entry'])

    # Filter gene and group entries
    genes = entries[entries['@type'] == 'gene'].copy()  # Filter gene entries
    genes['@name'] = genes['@name'].str.replace('hsa:', '', regex=False)  # Clean gene names

    groups = entries[entries['@type'] == 'group'].copy()  # Filter group entries

    # Create a mapping of gene IDs to names
    gene_id_name_map = genes.set_index('@id')['@name']  # Map IDs to names

    # Apply the mapping to extract gene names for each group
    groups['@name'] = groups.apply(lambda group: extract_group_genes(group, gene_id_name_map), axis=1)

    # Create a DataFrame for nodes combining genes and groups
    nodes = pd.DataFrame(columns=['id', 'type', 'genes'])  # Initialize nodes DataFrame
    gene_nodes = genes[['@id', '@type', '@name']].rename(columns={'@id': 'id', '@type': 'type', '@name': 'genes'})  # Rename columns for gene nodes
    group_nodes = groups[['@id', '@type', '@name']].rename(columns={'@id': 'id', '@type': 'type', '@name': 'genes'})  # Rename columns for group nodes
    nodes = pd.concat([gene_nodes, group_nodes], ignore_index=True)  # Concatenate gene and group nodes

    # Extract relations and reactions from the pathway
    relations = pd.DataFrame.from_dict(path_dict['pathway'].get('relation', []))  # Convert relations to DataFrame
    reactions = pd.DataFrame.from_dict(path_dict['pathway'].get('reaction', []))  # Convert reactions to DataFrame

    # Initialize the edges DataFrame to store relationships
    edges = pd.DataFrame(columns=['id1', 'id2', 'relation_type', 'relation_subtype'])  # Define columns for edges

    if not relations.empty:  # Check if there are any relations present
        # Filter relations for PPrel and GErel types
        pprel_rels = relations[(relations['@type'] == 'PPrel') | (relations['@type'] == 'GErel')]
        # Filter relations for PCrel type
        pcrel_rels = relations[relations['@type'] == 'PCrel']
    
        # Process PPrel and GErel relations
        for _, pprel in pprel_rels.iterrows():  # Iterate over each PPrel/GErel relation
            relation_type = pprel['@type']  # Get the type of relation
            relation_subtype = extract_relation_subtype(pprel.get('subtype'))  # Extract the subtype if available
            id1, id2 = pprel['@entry1'], pprel['@entry2']  # Get the IDs of the related entries
            # Create a new DataFrame for the edge
            new_edge = pd.DataFrame([[id1, id2, relation_type, relation_subtype]], columns=['id1', 'id2', 'relation_type', 'relation_subtype'])
            edges = pd.concat([edges, new_edge], ignore_index=True)  # Append the new edge to the edges DataFrame
    
        # Process PCrel relations
        for _, pcrel1 in pcrel_rels.iterrows():  # Iterate over each PCrel relation
            relation_type = pcrel1['@type']  # Get the type of relation
            relation_subtype = 'compound'  # Set the relation subtype to 'compound'
            if pcrel1['@entry1'] in nodes.values:  # Check if the first entry is a node
                id1 = pcrel1['@entry1']  # Get the ID of the first entry
                compound_id = pcrel1['@entry2']  # Get the ID of the compound entry
                # Retrieve all second entries that relate to the compound entry
                id2_list = pcrel_rels[pcrel_rels['@entry1'] == compound_id]['@entry2'].tolist()
                # Create a new DataFrame for multiple edges
                new_edges = pd.DataFrame([[id1, id2, relation_type, relation_subtype] for id2 in id2_list], columns=['id1', 'id2', 'relation_type', 'relation_subtype'])
                edges = pd.concat([edges, new_edges], ignore_index=True)  # Append new edges to the edges DataFrame

    if not reactions.empty:  # Check if there are any reactions present
        # Create a DataFrame to store gene-product relationships
        gene_product_substrate = pd.DataFrame(columns=['name', 'substrate', 'product', 'reaction_type'])
        
        # Iterate through the nodes DataFrame to access node details
        for node_index, node in nodes.iterrows():
            node_name = node['id']  # Access the node ID
            
            # Find the corresponding reaction in the reactions DataFrame
            reaction = reactions[reactions['@id'] == node_name]
            
            if not reaction.empty:  # Proceed if a reaction is found
                # Extract substrates and products from the matched reaction
                substrates = reaction['substrate'].iloc[0]  # List of substrate dictionaries
                products = reaction['product'].iloc[0]      # List of product dictionaries
                reaction_type = reaction['@type'].iloc[0]   # Extract reaction type
                
                # Process substrates and products into lists of IDs
                gene_product_substrate.loc[node_index, 'substrate'] = [s['@id'] for s in substrates] if isinstance(substrates, list) else [substrates['@id']]
                gene_product_substrate.loc[node_index, 'product'] = [p['@id'] for p in products] if isinstance(products, list) else [products['@id']]
                
                # Assign the gene name and reaction type
                gene_product_substrate.loc[node_index, 'name'] = node_name
                gene_product_substrate.loc[node_index, 'reaction_type'] = reaction_type  # Store the reaction type
        
        # Establish connections based on different associations
        gene_gene_list = []  # List to store gene-gene relationships
        
        for gene_index1 in gene_product_substrate.index:  # Iterate over each gene
            substrates1 = gene_product_substrate.loc[gene_index1, 'substrate']  # Get substrates of the first gene
            products1 = gene_product_substrate.loc[gene_index1, 'product']      # Get products of the first gene
            reaction_type1 = gene_product_substrate.loc[gene_index1, 'reaction_type']  # Get reaction type of the first gene
        
            # Iterate over products to find matches in other reactions
            for product in products1:
                for gene_index2 in gene_product_substrate.index:
                    if gene_index1 == gene_index2:  # Skip self-links
                        continue
        
                    substrates2 = gene_product_substrate.loc[gene_index2, 'substrate']  # Get substrates of the second gene
                    products2 = gene_product_substrate.loc[gene_index2, 'product']      # Get products of the second gene
                    reaction_type2 = gene_product_substrate.loc[gene_index2, 'reaction_type']  # Get reaction type of the second gene
        
                    # Condition 1: Check if product of reaction A is a substrate of reaction B
                    if product in substrates2:
                        final_relation_type = 'reversible' if reaction_type1 == 'reversible' and reaction_type2 == 'reversible' else 'irreversible'
                        gene_gene_list.append([
                            gene_product_substrate.loc[gene_index1, 'name'],
                            gene_product_substrate.loc[gene_index2, 'name'],
                            'substrate-product association',
                            final_relation_type
                        ])
        
                    # Condition 2: Check if product of reaction A is a product of reaction B (with reversible check)
                    elif product in products2 and reaction_type2 == 'reversible':
                        final_relation_type = 'reversible' if reaction_type1 == 'reversible' else 'irreversible'
                        gene_gene_list.append([
                            gene_product_substrate.loc[gene_index1, 'name'],
                            gene_product_substrate.loc[gene_index2, 'name'],
                            'product-product association',
                            final_relation_type
                        ])
        
            # Condition 3: If reaction A is reversible, check if substrates match with other reactions
            if reaction_type1 == 'reversible':
                for substrate in substrates1:
                    for gene_index2 in gene_product_substrate.index:
                        if gene_index1 == gene_index2:  # Skip self-links
                            continue
        
                        substrates2 = gene_product_substrate.loc[gene_index2, 'substrate']  # Get substrates of the second gene
                        reaction_type2 = gene_product_substrate.loc[gene_index2, 'reaction_type']  # Get reaction type of the second gene
        
                        # Check if substrates match
                        if substrate in substrates2:
                            final_relation_type = 'reversible' if reaction_type1 == 'reversible' and reaction_type2 == 'reversible' else 'irreversible'
                            gene_gene_list.append([
                                gene_product_substrate.loc[gene_index1, 'name'],
                                gene_product_substrate.loc[gene_index2, 'name'],
                                'substrate-substrate association',
                                final_relation_type
                            ])
        
        # Convert the gene-gene list to a DataFrame and append to edges
        gene_gene_df = pd.DataFrame(gene_gene_list, columns=['id1', 'id2', 'relation_type', 'relation_subtype'])
        
        # Assuming edges DataFrame is already defined and ready for appending
        edges = pd.concat([edges, gene_gene_df[['id1', 'id2', 'relation_type', 'relation_subtype']]], ignore_index=True)  # Append new edges to the edges DataFrame
    
    # Add bidirectional edges for specific subtypes
    bidirectional_subtypes = ['binding/association', 'reversible']
    for _, row in edges.iterrows():
        id1, id2, relation_subtype = row['id1'], row['id2'], row['relation_subtype']
        if relation_subtype in bidirectional_subtypes:
            reverse_edge = pd.DataFrame([[id2, id1, row['relation_type'], relation_subtype]], columns=['id1', 'id2', 'relation_type', 'relation_subtype'])
            edges = pd.concat([edges, reverse_edge], ignore_index=True)
    
    # Step 2: Create gene_set before grouping
    nodes['gene_set'] = nodes['genes'].apply(lambda x: frozenset(x.split(' ')))
    
    # Step 3: Group nodes by their genes and create a mapping for merged node IDs
    grouped_nodes = nodes.groupby('gene_set').agg({
        'id': 'first',  # Keep the first id for the group
        'type': 'first',  # Keep the first type for the group
        'genes': 'first'   # Keep the first gene name
    }).reset_index(drop=True)
    
    # Create a mapping of original IDs to the kept ID
    merged_id_mapping = {}
    for _, group in nodes.groupby('gene_set'):
        kept_id = group['id'].iloc[0]  # The ID to keep is the first one in the group
        for node_id in group['id']:
            merged_id_mapping[node_id] = kept_id  # Map each original ID to the kept ID
    
    # Update edges to replace merged node IDs
    edges['id1'] = edges['id1'].map(merged_id_mapping).fillna(edges['id1'])  # Replace id1 with mapped ID, retain original if no mapping
    edges['id2'] = edges['id2'].map(merged_id_mapping).fillna(edges['id2'])  # Replace id2 with mapped ID, retain original if no mapping
    
    return pathway_name, pathway_description, entries, grouped_nodes, edges

def create_network(nodes, edges):
    """
    Create a directed graph from nodes and edges DataFrames.

    Parameters:
    - nodes (DataFrame): Contains node IDs.
    - edges (DataFrame): Contains directed edges with 'id1' and 'id2'.

    Returns:
    - network (DiGraph): Directed graph of nodes and edges.
    """
    
    network = nx.DiGraph()  # Initialize a directed graph
    network.add_nodes_from(nodes['id'])  # Add nodes from the 'id' column
    network.add_edges_from([(row['id1'], row['id2']) for _, row in edges.iterrows()])  # Add edges as tuples (id1, id2)
    
    return network

def plot_network(network, nodes, pathway_name, pathway_description, plot_dir=None):
    """
    Visualize a directed network and save the plot as a PNG file.

    Parameters:
    - network (DiGraph): The directed graph representing the network.
    - nodes (DataFrame): DataFrame containing node attributes, including type and gene names.
    - pathway_name (str): The name of the pathway.
    - pathway_description (str): A description of the pathway.
    - plot_dir (str, optional): Directory to save the plot; if None, the plot is not saved.
    """
    plt.figure(figsize=(16, 12))  # Set figure size
    pos = nx.spring_layout(network, k=0.6)  # Generate layout for nodes

    node_colors = []  # List to store colors for each node
    node_labels = {}  # Dictionary to store labels for each node

    for node in network.nodes:  # Iterate through each node in the network
        node_data = nodes[nodes['id'] == node]  # Get node data from DataFrame
        if not node_data.empty:  # Check if node data is available
            node_type = node_data['type'].values[0]  # Get the type of the node
            gene_names = node_data['genes'].values[0]  # Get gene names associated with the node
        else:
            node_type = 'unknown'  # Default type if no data
            gene_names = 'Unknown'  # Default gene name

        # Assign colors based on node type
        node_colors.append('lightblue' if node_type == 'gene' else 'lightcoral' if node_type == 'group' else 'lightblue')
        node_labels[node] = gene_names.split(' ')[0]  # Use the first gene name as the label

    # Draw network components
    nx.draw_networkx_nodes(network, pos, node_size=600, node_color=node_colors, alpha=0.7)
    nx.draw_networkx_edges(network, pos, width=3, alpha=0.8, edge_color='gray', arrowstyle='-|>')
    nx.draw_networkx_labels(network, pos, labels=node_labels, font_size=11, font_color='black')

    plt.title(f"Network: {pathway_name} - {pathway_description}")  # Set plot title
    plt.axis('off')  # Hide axes

    # Save plot if directory is specified
    if plot_dir:
        filename = os.path.join(plot_dir, f"{pathway_name}_network.png")
        plt.savefig(filename, bbox_inches='tight', dpi=300)  # Save as PNG
    plt.close()  # Close the plot

def get_node_centrality(network, pathway_name, pathway_description):
    """
    Calculate centrality metrics for nodes in a directed network.

    Parameters:
    - network (DiGraph): The directed graph representing the network.
    - pathway_name (str): The name of the pathway.
    - pathway_description (str): A description of the pathway.

    Returns:
    - DataFrame: A DataFrame containing centrality metrics for each node, including degree, closeness, and betweenness.
    """
    # Check if the network is empty
    if len(network.nodes()) == 0:
        return pd.DataFrame({
            'Pathway_Name': pathway_name,
            'Description': pathway_description,
            'id': [],
            'Degree': [],
            'In_Degree': [],
            'Out_Degree': [],
            'Closeness': [],
            'Local_Reaching': [],
            'Betweenness': []
        })

    # Check if the network has no edges
    if len(network.edges()) == 0:
        num_nodes = len(network.nodes())
        return pd.DataFrame({
            'Pathway_Name': [pathway_name] * num_nodes,
            'Description': [pathway_description] * num_nodes,
            'id': list(network.nodes()),
            'Degree': [0] * num_nodes,
            'In_Degree': [0] * num_nodes,
            'Out_Degree': [0] * num_nodes,
            'Closeness': [0] * num_nodes,
            'Local_Reaching': [0] * num_nodes,
            'Betweenness': [0] * num_nodes
        })

    # Initialize centrality dictionaries
    degree_dict = dict(network.degree())  # Calculate degree for each node
    betweenness_dict = nx.betweenness_centrality(network)  # Calculate betweenness centrality
    
    # Compute in-degree and out-degree (only for directed graphs)
    if nx.is_directed(network):
        in_degree_dict = dict(network.in_degree())  # Calculate in-degree for each node
        out_degree_dict = dict(network.out_degree())  # Calculate out-degree for each node
    else:
        in_degree_dict = degree_dict  # For undirected graphs, in-degree equals degree
        out_degree_dict = degree_dict  # For undirected graphs, out-degree equals degree
    
    # Compute closeness centrality
    closeness_dict = nx.closeness_centrality(network)  # Calculate closeness centrality
    
    # Compute local reaching centrality for all nodes
    local_reaching_dict = {}
    for node in network.nodes():
        local_reaching_dict[node] = nx.local_reaching_centrality(network, v=node)

    # Create DataFrame for node centralities
    node_centrality_df = pd.DataFrame({
        'Pathway_Name': pathway_name,
        'Description': pathway_description,
        'id': list(degree_dict.keys()),
        'Degree': list(degree_dict.values()),
        'In_Degree': [in_degree_dict.get(node, 0) for node in degree_dict.keys()],
        'Out_Degree': [out_degree_dict.get(node, 0) for node in degree_dict.keys()],
        'Closeness': [closeness_dict.get(node, 0) for node in degree_dict.keys()],
        'Local_Reaching': [local_reaching_dict.get(node, 0) for node in degree_dict.keys()],
        'Betweenness': [betweenness_dict.get(node, 0) for node in degree_dict.keys()]  # Default to 0 for isolated nodes
    })

    return node_centrality_df  

def get_gene_centrality(node_centrality_df, nodes):
    # Map node IDs to corresponding genes
    id_gene_map = nodes.set_index('id')['genes']
    
    # Merge node centrality data with gene mappings
    node_centrality_df['genes'] = node_centrality_df['id'].map(id_gene_map)
    
    # Explode 'genes' column to give each gene the same centrality value
    exploded_df = node_centrality_df.assign(genes=node_centrality_df['genes'].str.split(' ')).explode('genes')
    
    # Group by gene, sum centrality values, and retain Pathway_Name and Description
    gene_centrality_df = exploded_df.groupby('genes').agg({
        'Pathway_Name': 'first',      # Retain pathway name
        'Description': 'first',       # Retain pathway description
        'Degree': 'sum',              # Sum of degree centrality
        'In_Degree': 'sum',           # Sum of in-degree
        'Out_Degree': 'sum',          # Sum of out-degree
        'Closeness': 'sum',           # Sum of closeness centrality
        'Local_Reaching': 'sum',      # Sum of local reaching centrality
        'Betweenness': 'sum'          # Sum of betweenness centrality
    }).reset_index()

    return gene_centrality_df

def kegg_topology(kgml_file, output_dir, plot):
    """
    Analyze the KEGG pathway from a KGML file, compute node centrality metrics, and optionally plot the network.

    Parameters:
    - kgml_file (str): Path to the KGML file containing pathway data.
    - output_dir (str): Directory to save output files and plots.
    - plot (bool): Whether to generate and save a plot of the pathway network.

    Returns:
    - DataFrame: Centrality metrics for genes in the pathway, including Degree, In_Degree, Out_Degree, Closeness, Local_Reaching, and Betweenness.
    - nodes: Node information extracted from the pathway.
    - edges: Edge information extracted from the pathway.
    """
    # Open the KGML file and parse its contents into a dictionary
    with open(kgml_file, 'r') as file:
        path_dict = xmltodict.parse(file.read())
    
    # Extract pathway name, description, nodes, and edges from the parsed dictionary
    pathway_name, pathway_description, entries, nodes, edges = get_nodes_and_edges(path_dict)
    
    # Create a directed network from the nodes and edges
    network = create_network(nodes, edges)
    
    # Check if the network is empty; if so, skip further processing
    if len(network.nodes()) == 0:
        print(f"Skipping empty graph for {kgml_file}.")
        return pd.DataFrame(columns=['Pathway_Name', 'Description', 'genes', 'Degree', 'In_Degree', 'Out_Degree', 'Closeness', 'Local_Reaching', 'Betweenness'])
    
    # Compute centrality metrics for the nodes in the network
    node_centrality_df = get_node_centrality(network, pathway_name, pathway_description)
    
    # Compute gene centrality metrics based on the node centrality DataFrame
    gene_centrality_df = get_gene_centrality(node_centrality_df, nodes)
    
    # If plotting is enabled, create and save the network plot
    if plot:
        # Create a directory for plots within the specified output directory
        plot_dir = os.path.join(output_dir, 'pathway_graphs/')
        os.makedirs(plot_dir, exist_ok=True)  # Create directory if it doesn't exist
        
        # Plot the network and save the figure
        plot_network(network, nodes, pathway_name, pathway_description, plot_dir)
    
    return gene_centrality_df, nodes, edges  # Return the centrality DataFrame and node/edge information


def load_mapping(mapping_file):
    # Load the mapping file, ensuring it's in the correct format
    mapping = pd.read_csv(mapping_file, sep='\t', header=None, dtype=str)

    # Check if the mapping file has exactly two columns
    if mapping.shape[1] != 2:
        raise ValueError(f"Expected 2 columns, but found {mapping.shape[1]}.")

    # Rename columns to a consistent lowercase for case-insensitive checks
    mapping.columns = ['ID_1', 'ID_2']
    id_1_lower = mapping['ID_1'].str.lower().iloc[0]
    id_2_lower = mapping['ID_2'].str.lower().iloc[0]

    # Rename columns based on content
    if 'ensembl' in id_2_lower:
        mapping.columns = ['ENTREZID', 'ENSEMBL']
    elif 'entrez' in id_1_lower:
        mapping.columns = ['ENSEMBL', 'ENTREZID']

    # Clean up: remove NaNs and duplicates, strip whitespace
    mapping = (mapping.dropna()
                     .drop_duplicates()
                     .assign(ENTREZID=lambda df: df['ENTREZID'].str.strip(),
                             ENSEMBL=lambda df: df['ENSEMBL'].str.strip()))

    # Return a dictionary mapping ENTREZID to ENSEMBL
    return dict(zip(mapping['ENTREZID'], mapping['ENSEMBL']))

def replace_entrez_with_ensembl(df, mapping_dict):
    # Replace Entrez IDs in 'Entrez_ID' column with Ensembl IDs
    df['Ensembl_ID'] = df['Entrez_ID'].apply(
        lambda x: ', '.join([mapping_dict.get(gene, '') for gene in x.split(', ')])
    )
    # Remove rows where Ensembl_ID is empty
    df = df[df['Ensembl_ID'].str.strip() != '']
    return df

def kegg_topology_all(kgml_dir, output_dir, mapping_file=None, plot=False):
    os.makedirs(output_dir, exist_ok=True)
    
    kegg_centrality_df = pd.DataFrame(columns=['Pathway_Name', 'Description', 'genes', 
                                                'Degree', 'In_Degree', 'Out_Degree', 
                                                'Closeness', 'Local_Reaching', 'Betweenness'])
    
    kgml_files = [os.path.join(kgml_dir, file) for file in os.listdir(kgml_dir) if file.endswith('.kgml')]
    
    def process_file(kgml_file):
        gene_centrality_df, _, _ = kegg_topology(kgml_file, output_dir, plot)
        return gene_centrality_df

    results = Parallel(n_jobs=-1)(
        delayed(process_file)(kgml_file) for kgml_file in tqdm(kgml_files, desc="Processing KGML files")
    )

    kegg_centrality_df = pd.concat(results, ignore_index=True)

    # Save the DataFrame with Entrez IDs (renaming genes to Entrez_ID)
    kegg_centrality_df.rename(columns={'genes': 'Entrez_ID'}, inplace=True)
    output_file_entrez = os.path.join(output_dir, 'kegg_centrality_results_entrez.tsv')
    kegg_centrality_df[['Pathway_Name', 'Description', 'Entrez_ID', 
                             'Degree', 'In_Degree', 'Out_Degree', 
                             'Closeness', 'Local_Reaching', 'Betweenness']].to_csv(output_file_entrez, index=False,sep = '\t')
    
    if mapping_file:
        mapping_dict = load_mapping(mapping_file)

        # Create a DataFrame with Ensembl IDs and drop the Entrez_ID column
        kegg_centrality_df = replace_entrez_with_ensembl(kegg_centrality_df, mapping_dict)
        output_file_ensembl = os.path.join(output_dir, 'kegg_centrality_results_ensembl.tsv')
        kegg_centrality_df[['Pathway_Name', 'Description', 'Ensembl_ID', 
                             'Degree', 'In_Degree', 'Out_Degree', 
                             'Closeness', 'Local_Reaching', 'Betweenness']].to_csv(output_file_ensembl, index=False,sep = '\t')

    
    return kegg_centrality_df


def parse_kgml_main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Process KEGG pathway topology.")
    
    # Add arguments
    parser.add_argument('--kgml_dir', type=str, required=True, 
                        help='Directory containing KGML files.')
    parser.add_argument('--output_dir', type=str, required=True, 
                        help='Directory to save output files.')
    parser.add_argument('--mapping_file', type=str, 
                        help='File mapping ENTREZ IDs to ENSEMBL IDs.')
    parser.add_argument('--plot', action='store_true', 
                        help='Flag to enable plotting.')
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # Call the kegg_topology_all function
    kegg_topology_all(args.kgml_dir, args.output_dir, 
                                       mapping_file=args.mapping_file, 
                                       plot=args.plot)


if __name__ == '__main__':
    #parse_kgml_main()
    # Uncomment the previous line and comment the following lines to run from the command line
    
    # Direct invocation for IDE usage
    output_dir = '../../../data/pathways/KEGG/topology/'
    kgml_dir = '../../../data/pathways/KEGG/kgml'
    mapping_file = '../../../data/id_conversion/ENTREZID_to_ENSEMBL.txt'
    centrality_df = kegg_topology_all(kgml_dir, output_dir, mapping_file=mapping_file, plot=True)

