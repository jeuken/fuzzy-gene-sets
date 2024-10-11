#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
from Bio.KEGG.KGML.KGML_parser import read
import networkx as nx
import matplotlib.pyplot as plt

# Define the file paths
kgml_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/kgml/hsa00010.kgml"
reactions_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/reactions.csv"
relations_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/relations.csv"

# Load the pathway from the KGML file
with open(kgml_path, 'r') as kgml_file:
    pathway = read(kgml_file)

# Save reactions to a CSV file
with open(reactions_path, 'w', newline='') as reactions_file:
    writer = csv.writer(reactions_file)
    writer.writerow(["Reaction ID", "Reaction Type", "Substrates", "Products"])
    for reaction in pathway.reactions:
        substrates = "; ".join([s.name for s in reaction.substrates])
        products = "; ".join([p.name for p in reaction.products])
        writer.writerow([reaction.id, reaction.type, substrates, products])

print(f"Reactions saved to {reactions_path}")

# Save relations to a CSV file
with open(relations_path, 'w', newline='') as relations_file:
    writer = csv.writer(relations_file)
    writer.writerow(["Entry1", "Entry2", "Relation Type", "Subtype"])
    for relation in pathway.relations:
        subtypes = "; ".join([f"{subtype[0]} ({subtype[1]})" for subtype in relation.subtypes])
        writer.writerow([relation.entry1.name, relation.entry2.name, relation.type, subtypes])

print(f"Relations saved to {relations_path}")

# Initialize a NetworkX directed graph
G = nx.DiGraph()

# Add nodes and edges for reactions
for reaction in pathway.reactions:
    for substrate in reaction.substrates:
        G.add_node(substrate.name, type='substrate')
    for product in reaction.products:
        G.add_node(product.name, type='product')
        for substrate in reaction.substrates:
            G.add_edge(substrate.name, product.name, reaction_type=reaction.type)

# Add nodes and edges for relations
for relation in pathway.relations:
    entry1 = relation.entry1.name
    entry2 = relation.entry2.name
    G.add_node(entry1, type='relation')
    G.add_node(entry2, type='relation')
    for subtype in relation.subtypes:
        G.add_edge(entry1, entry2, relation_type=relation.type, subtype=subtype[0])

# Visualize the graph
plt.figure(figsize=(12, 8))
pos = nx.spring_layout(G)  # Layout for better visualization
nx.draw(G, pos, with_labels=True, node_size=500, node_color="lightblue", 
        font_size=8, font_weight="bold", edge_color="gray")
plt.title("KEGG Pathway Visualization")
plt.show()

