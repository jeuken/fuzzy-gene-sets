import os
import xmltodict
import pandas as pd
import networkx as nx  
import matplotlib.pyplot as plt  
from joblib import Parallel, delayed  
from tqdm import tqdm

# Extracts the '@name' from a subtype input, which can be a list or a dictionary.
def extract_subtype(subtype):
    # If subtype is a non-empty list, return '@name' of the first element if it exists
    if isinstance(subtype, list) and subtype:
        return subtype[0].get('@name')
    # If subtype is a dictionary, return its '@name' if it exists
    elif isinstance(subtype, dict):
        return subtype.get('@name')
    # Return subtype as is if it doesn't match the above conditions
    return subtype


# Function to process each KGML file
def process_kgml_file(kgml_file):
    # Open the specified KGML file and parse its content into a dictionary
    with open(kgml_file, 'r') as file:
        path_dict = xmltodict.parse(file.read())

    # Extract pathway name and description, cleaning the name
    pathway_name = path_dict['pathway']['@name'].replace('path:', '')
    pathway_description = path_dict['pathway'].get('@title', '')

    # Convert entries to a DataFrame and filter for genes
    entries = pd.DataFrame.from_dict(path_dict['pathway']['entry'])  
    genes = entries[entries['@type'] == 'gene']
    
    # Clean gene names by removing the 'hsa:' prefix
    genes.loc[:, '@name'] = genes['@name'].str.replace('hsa:', '', regex=False)  
    
    # Create DataFrames for relations and reactions, if they exist
    relations = pd.DataFrame.from_dict(path_dict['pathway'].get('relation', []))  
    reactions = pd.DataFrame.from_dict(path_dict['pathway'].get('reaction', [])) 
  
    # Initialize an empty DataFrame to store gene-gene relationships
    gene_gene = pd.DataFrame(columns=['gene1', 'gene2', 'id1','id2','subtype'])
    
    # Create a mapping from gene IDs to gene names for easy reference
    id_to_name = dict(zip(genes['@id'], genes['@name']))  

    # Filter relations to include only gene pairs
    if not relations.empty:
        gene_ids = set(genes['@id'])
        
        pprel_rels = relations[(relations['@type'] == 'PPrel') | (relations['@type'] == 'GErel')]
        ecrel_rels = relations[relations['@type'] == 'ECrel']
        pcrel_rels = relations[relations['@type'] == 'PCrel']
        
        # Process PPrel relations
        for _, pprel in pprel_rels.iterrows():
            entry1, entry2 = pprel['@entry1'], pprel['@entry2']
        
            if entry1 in gene_ids and entry2 in gene_ids:
                gene1_name = id_to_name.get(entry1, entry1)
                gene2_name = id_to_name.get(entry2, entry2)
                subtype = extract_subtype(pprel.get('subtype'))
                gene_gene = pd.concat([gene_gene, pd.DataFrame([[gene1_name, gene2_name,entry1,entry2, subtype]],
                                                               columns=['gene1', 'gene2','id1','id2', 'subtype'])], ignore_index=True)
        
            # Handle gene-group relations (gene to group)
            elif entry1 in gene_ids and entry2 in entries[entries['@type'] == 'group']['@id'].values:
                gene1_name = id_to_name.get(entry1, entry1)
                group_genes = entries.loc[entries['@id'] == entry2, 'component'].values[0]
                group_gene_ids = [gene['@id'] for gene in group_genes if isinstance(gene, dict) and '@id' in gene]
                group_gene_names = [id_to_name.get(group_gene_id, group_gene_id) for group_gene_id in group_gene_ids]
                group_gene_name_str = ' '.join(group_gene_names)
                subtype = extract_subtype(pprel.get('subtype'))
                
                gene_gene = pd.concat([gene_gene, pd.DataFrame([[gene1_name, group_gene_name_str,entry1,entry2, subtype]],
                                                                columns=['gene1', 'gene2','id1','id2', 'subtype'])], ignore_index=True)
            # Handle group-gene relations (group to gene)
            elif entry1 in entries[entries['@type'] == 'group']['@id'].values and entry2 in gene_ids:
                gene2_name = id_to_name.get(entry2, entry2)
                group_genes = entries.loc[entries['@id'] == entry1, 'component'].values[0]
                group_gene_ids = [gene['@id'] for gene in group_genes if isinstance(gene, dict) and '@id' in gene]
                group_gene_names = [id_to_name.get(group_gene_id, group_gene_id) for group_gene_id in group_gene_ids]
                group_gene_name_str = ' '.join(group_gene_names)
        
                subtype = extract_subtype(pprel.get('subtype'))
                
                gene_gene = pd.concat([gene_gene, pd.DataFrame([[group_gene_name_str, gene2_name,entry1,entry2, subtype]],
                                                                columns=['gene1', 'gene2','id1','id2', 'subtype'])], ignore_index=True)
        
            # Handle group-group relations (group to group)
            elif entry1 in entries[entries['@type'] == 'group']['@id'].values and entry2 in entries[entries['@type'] == 'group']['@id'].values:
                group1_genes = entries.loc[entries['@id'] == entry1, 'component'].values[0]
                group2_genes = entries.loc[entries['@id'] == entry2, 'component'].values[0]
        
                group1_gene_ids = [gene['@id'] for gene in group1_genes if isinstance(gene, dict) and '@id' in gene]
                group2_gene_ids = [gene['@id'] for gene in group2_genes if isinstance(gene, dict) and '@id' in gene]
                
                group1_gene_names = [id_to_name.get(gene_id, gene_id) for gene_id in group1_gene_ids]
                group2_gene_names = [id_to_name.get(gene_id, gene_id) for gene_id in group2_gene_ids]
        
                group1_gene_name_str = ' '.join(group1_gene_names)
                group2_gene_name_str = ' '.join(group2_gene_names)
        
                subtype = extract_subtype(pprel.get('subtype'))
                
                gene_gene = pd.concat([gene_gene, pd.DataFrame([[group1_gene_name_str, group2_gene_name_str,entry1,entry2, subtype]],
                                                                columns=['gene1', 'gene2','id1','id2', 'subtype'])], ignore_index=True)


        # Process entries
        for _, ecrel in ecrel_rels.iterrows():
            entry1 = ecrel['@entry1']
            entry2 = ecrel['@entry2']
        
            # Check if '@reaction' exists before accessing it
            if '@reaction' in entries.columns:
                reaction1 = entries[entries['@id'] == entry1]['@reaction'].values
                reaction2 = entries[entries['@id'] == entry2]['@reaction'].values
        
                # Get the corresponding subtype for each reaction
                if reaction1.size > 0:
                    subtype1 = reactions[reactions['@name'] == reaction1[0]]['@type'].values
                else:
                    subtype1 = []
        
                if reaction2.size > 0:
                    subtype2 = reactions[reactions['@name'] == reaction2[0]]['@type'].values
                else:
                    subtype2 = []
        
                # Determine subtype
                if subtype1.size > 0 and subtype2.size > 0:
                    if (subtype1[0] == "reversible") and (subtype2[0] == "reversible"):
                        subtype = "reversible"
                    else:
                        subtype = "irreversible"
                else:
                    subtype = "unknown"  # Handle missing subtypes if necessary
            else:
                subtype = "unknown" 
                print("Column '@reaction' not found in entries DataFrame")
                # Handle this case as needed
            
            # Check if entry1 is a gene and entry2 is a gene
            if entry1 in genes['@id'].values and entry2 in genes['@id'].values:
                gene1_name = id_to_name.get(entry1, entry1)
                gene2_name = id_to_name.get(entry2, entry2)
                
                # Check for subtype and extract name if available
                subtype = subtype
                
                gene_gene = pd.concat([gene_gene, pd.DataFrame([[gene1_name, gene2_name, entry1,entry2, subtype]],
                                                               columns=['gene1', 'gene2','id1','id2', 'subtype'])], ignore_index=True)
        
            elif entry1 in genes['@id'].values and entry2 in entries[entries['@type'] == 'group']['@id'].values:
                gene1_name = id_to_name.get(entry1, entry1)
                group_genes = entries.loc[entries['@id'] == entry2, 'component'].values[0]
                group_gene_ids = [gene['@id'] for gene in group_genes if isinstance(gene, dict) and '@id' in gene]
                group_gene_names = [id_to_name.get(group_gene_id, group_gene_id) for group_gene_id in group_gene_ids]
                group_gene_name_str = ' '.join(group_gene_names)
                
                gene_gene = pd.concat([gene_gene, pd.DataFrame([[gene1_name, group_gene_name_str,entry1,entry2, subtype]],
                                                                columns=['gene1', 'gene2','id1','id2', 'subtype'])], ignore_index=True)
            # Check if entry1 is a group and entry2 is a gene
            elif entry1 in entries[entries['@type'] == 'group']['@id'].values and entry2 in genes['@id'].values:
                gene2_name = id_to_name.get(entry2, entry2)
                group_genes = entries.loc[entries['@id'] == entry1, 'component'].values[0]
                group_gene_ids = [gene['@id'] for gene in group_genes if isinstance(gene, dict) and '@id' in gene]
                group_gene_names = [id_to_name.get(group_gene_id, group_gene_id) for group_gene_id in group_gene_ids]
                group_gene_name_str = ' '.join(group_gene_names)
                
                gene_gene = pd.concat([gene_gene, pd.DataFrame([[group_gene_name_str, gene2_name,entry1,entry2, subtype]],
                                                                columns=['gene1', 'gene2','id1','id2', 'subtype'])], ignore_index=True)
        

            # Check if both entry1 and entry2 are groups
            elif entry1 in entries[entries['@type'] == 'group']['@id'].values and entry2 in entries[entries['@type'] == 'group']['@id'].values:
                group1_genes = entries.loc[entries['@id'] == entry1, 'component'].values[0]
                group2_genes = entries.loc[entries['@id'] == entry2, 'component'].values[0]
        
                group1_gene_ids = [gene['@id'] for gene in group1_genes if isinstance(gene, dict) and '@id' in gene]
                group2_gene_ids = [gene['@id'] for gene in group2_genes if isinstance(gene, dict) and '@id' in gene]
                
                group1_gene_names = [id_to_name.get(gene_id, gene_id) for gene_id in group1_gene_ids]
                group2_gene_names = [id_to_name.get(gene_id, gene_id) for gene_id in group2_gene_ids]
        
                group1_gene_name_str = ' '.join(group1_gene_names)
                group2_gene_name_str = ' '.join(group2_gene_names)
                
                gene_gene = pd.concat([gene_gene, pd.DataFrame([[group1_gene_name_str, group2_gene_name_str,entry1,entry2, subtype]],
                                                                columns=['gene1', 'gene2','id1','id2', 'subtype'])], ignore_index=True)

        # Handle PCrel to find indirect relations
          # Select PCrel relations
        for _, pcrel1 in pcrel_rels.iterrows():
            if pcrel1['@entry1'] in gene_ids:
                gene1 = pcrel1['@entry1']
                compound_id = pcrel1['@entry2']
                associated_genes = pcrel_rels[pcrel_rels['@entry1'] == compound_id]['@entry2'].tolist()
                indirect_rels = [
                    [id_to_name.get(gene1, gene1), id_to_name.get(gene2, gene2), gene1,gene2,'indirect_relation']
                    for gene2 in associated_genes
                ]
                gene_gene = pd.concat([gene_gene, pd.DataFrame(indirect_rels, columns=['gene1', 'gene2','id1','id2', 'subtype'])], ignore_index=True)

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
                    group_gene_list.append([gene1_name, gene2_name,gene_ids[i],gene_ids[j], 'group_association'])

        if group_gene_list:  # Create DataFrame if there are group relations
            gene_gene = pd.concat([gene_gene, pd.DataFrame(group_gene_list, columns=['gene1', 'gene2','id1','id2', 'subtype'])], ignore_index=True)


     # Add bidirectional relations for specified subtypes
    bidirectional_subtypes = ['binding/association', 'group_association','reversible']
    for _, row in gene_gene.iterrows():
        gene1 = row['gene1']
        gene2 = row['gene2']
        id1 = row['id1']
        id2 = row['id2']
        subtype = row['subtype']
         
        # If subtype is bidirectional, add the reverse relation
        if subtype in bidirectional_subtypes:
            gene_gene = pd.concat([gene_gene, pd.DataFrame([[gene2, gene1, id2,id1, subtype]], columns=['gene1', 'gene2','id1','id2', 'subtype'])], ignore_index=True)

    print(gene_gene)
    #gene_gene['gene1'] = gene_gene['gene1'].apply(lambda x: x.split() if isinstance(x, str) else [x])
    #gene_gene['gene2'] = gene_gene['gene2'].apply(lambda x: x.split() if isinstance(x, str) else x)
    
    # Explode the gene_gene DataFrame to have individual gene entries
    #gene_gene_explode = gene_gene.explode('gene1').explode('gene2').reset_index(drop=True)  # Explode both gene1 and gene2

    # Create directed graph
    G = nx.DiGraph()
    
    edges = tuple(map(tuple, gene_gene[['gene1', 'gene2']].values))
    
    # Add edges from the gene_gene DataFrame
    G.add_edges_from(edges)
    
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
        'Betweenness_Membership': [betweenness_dict.get(gene, 0) for gene in degree_dict.keys()],  # Default to 0 for isolated nodes
        'Pathway_Name': pathway_name,  # Add pathway name
        'Description': pathway_description  # Add pathway description
    })
    
    # Explode gene names if they contain multiple IDs
    centrality_df['gene'] = centrality_df['gene'].apply(lambda x: x.split())
    
    # Explode the 'gene' column to create a long format table where each gene has its own row
    long_centrality_df = centrality_df.explode('gene').reset_index(drop=True)
    
    # Group by 'gene', sum the relevant metrics (degree and Betweenness_Membership), and take the first value for 'Pathway_Name' and 'Description'
    long_centrality_df = long_centrality_df.groupby('gene', as_index=False).agg({
        'degree': 'sum',
        'Betweenness_Membership': 'sum',
        'Pathway_Name': 'first',  # Take the first occurrence (assuming it's the same for all rows)
        'Description': 'first'    # Take the first occurrence (assuming it's the same for all rows)
    })
    
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
    plt.title(f"Gene Interaction Network from {pathway_name}: {pathway_description}")
    
    plt.axis('off')
    
    # Save the plot with the pathway name and ID
    plt.savefig(os.path.join(plot_directory, file_name), format='png',dpi=300)
    plt.close()  
    
    return long_centrality_df

# Use joblib for parallel processing
pathway_files = [os.path.join('../../data/pathways/KEGG/kgml/', file) for file in os.listdir('../../data/pathways/KEGG/kgml/') if file.endswith('.kgml')]
results = Parallel(n_jobs=-1)(delayed(process_kgml_file)(file) for file in tqdm(pathway_files))

# Concatenate results into a single DataFrame and handle empty DataFrames
final_results = pd.concat(results, ignore_index=True).dropna().reset_index(drop=True)

# Save the final results DataFrame to a CSV file
final_results.to_csv('../../data/pathways/KEGG/Centrality_Results.csv', index=False)