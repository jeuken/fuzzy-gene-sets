import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import os
import json

# Define the input file path
input_file = '/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/kegg_gene_network_hsa/hsa00020.tsv'

# Load the edges from the file
edges = pd.read_csv(input_file, sep='\t')

# Create a directed graph from the edge list
graph = nx.from_pandas_edgelist(edges, source='entry1', target='entry2', create_using=nx.DiGraph())

# Load the node positions from the JSON file
#with open('/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/kegg_gene_network_hsa/hsa05414_graphics.txt') as file:
#    pos = json.load(file)

# Plot the graph with node labels
plt.figure(figsize=(12, 8))

# Customize the appearance
nx.draw_networkx(graph, with_labels=True, node_size=500, node_color='lightblue', 
                 font_size=8, font_color='black', font_weight='bold', edge_color='gray')

# Adjust layout
plt.tight_layout()

# Extract the pathway name from the input file (without extension)
pathway_name = os.path.basename(input_file).split('.')[0]

# Define the output directory and file path
output_dir = '/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/pathway_plots/'
output_file = os.path.join(output_dir, f'{pathway_name}_knext.png')

# Save the figure with the pathway name and high resolution
plt.savefig(output_file, format='png', dpi=300)

# Show the plot
plt.show()
