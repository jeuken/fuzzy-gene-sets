#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
from Bio.KEGG.KGML.KGML_parser import read

# Load the pathway from the KGML file
file_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/kgml/hsa00010.kgml"
with open(file_path, 'r') as kgml_file:
    pathway = read(kgml_file)

# Save reactions to a CSV file
reactions_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/reactions.csv"
with open(reactions_path, 'w', newline='') as reactions_file:
    writer = csv.writer(reactions_file)
    writer.writerow(["Reaction ID", "Reaction Type", "Substrates", "Products"])
    
    for reaction in pathway.reactions:
        substrates = "; ".join([s.name for s in reaction.substrates])
        products = "; ".join([p.name for p in reaction.products])
        writer.writerow([reaction.id, reaction.type, substrates, products])

print(f"Reactions saved to {reactions_path}")

# Save relations to a CSV file
relations_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/relations.csv"
with open(relations_path, 'w', newline='') as relations_file:
    writer = csv.writer(relations_file)
    writer.writerow(["Entry1", "Entry2", "Relation Type", "Subtype"])
    
    for relation in pathway.relations:
        subtypes = "; ".join([f"{subtype[0]} ({subtype[1]})" for subtype in relation.subtypes])
        writer.writerow([relation.entry1.name, relation.entry2.name, relation.type, subtypes])

print(f"Relations saved to {relations_path}")
