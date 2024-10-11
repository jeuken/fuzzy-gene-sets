#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from Bio.KEGG.KGML.KGML_parser import read
from Bio.KEGG.KGML.KGML_pathway import Pathway, Entry, Reaction, Relation

def load_and_parse_kgml_files(directory):
    """
    Load and parse KGML files from the specified directory.

    Args:
        directory (str): The directory containing the KGML files.

    Returns:
        dict: A dictionary with pathway IDs as keys and extracted information as values.
    """
    pathway_data = {}

    # Iterate over all KGML files in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".kgml"):
            pathway_id = filename[:-5]  # Remove '.kgml' extension
            file_path = os.path.join(directory, filename)
            
            # Parse the KGML file
            with open(file_path, "r") as file:
                pathway = read(file)
            
            # Extract relevant information
            entries = pathway.entries
            reactions = pathway.reactions
            relations = pathway.relations
            
            # Prepare data for the pathway
            entry_data = [{"ID": entry.id, "Type": entry.type, "Names": entry.names} for entry in entries.values()]
            reaction_data = [{"ID": reaction.id, "Type": reaction.type,
                              "Substrates": [sub.id for sub in reaction.substrates],
                              "Products": [prod.id for prod in reaction.products]} for reaction in reactions]
            relation_data = [{"Entry1": relation.entry1, "Entry2": relation.entry2, "Type": relation.type}
                             for relation in relations]
            
            # Store the data
            pathway_data[pathway_id] = {
                "entries": entry_data,
                "reactions": reaction_data,
                "relations": relation_data
            }

    return pathway_data

# Example usage
if __name__ == "__main__":
    directory = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/kgml"
    data = load_and_parse_kgml_files(directory)
    
    # Print out the first pathway's data as an example
    first_pathway_id = next(iter(data))
    print(f"Data for pathway {first_pathway_id}:")
    print(data[first_pathway_id])


