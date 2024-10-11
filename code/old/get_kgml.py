#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import requests

# KEGG REST API URL for pathways
KEGG_BASE_URL = "https://rest.kegg.jp"

# Directory to save downloaded KGML files
OUTPUT_DIR = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/kgml"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Function to get a list of all pathway IDs
def get_all_pathways(organism="hsa"):
    response = requests.get(f"{KEGG_BASE_URL}/list/pathway/{organism}")
    if response.status_code == 200:
        pathways = response.text.splitlines()
        pathway_ids = [line.split('\t')[0] for line in pathways if '\t' in line]  # Extract pathway IDs
        return pathway_ids
    else:
        print("Failed to retrieve pathway list")
        return []

# Function to download KGML file for a specific pathway
def download_kgml(pathway_id):
    url = f"{KEGG_BASE_URL}/get/{pathway_id}/kgml"
    response = requests.get(url)
    if response.status_code == 200:
        file_path = os.path.join(OUTPUT_DIR, f"{pathway_id}.kgml")
        with open(file_path, "w") as file:
            file.write(response.text)
        print(f"Downloaded and saved KGML file for pathway: {pathway_id}")
    else:
        print(f"Failed to download KGML for pathway {pathway_id}")

# Main function to download KGML files
def main():
    # Get all pathway IDs for human (hsa)
    pathway_ids = get_all_pathways()

    # Download each KGML file
    for index, pathway_id in enumerate(pathway_ids):
        print(f"Downloading KGML file for pathway {index + 1} of {len(pathway_ids)}: {pathway_id}")
        download_kgml(pathway_id)

if __name__ == "__main__":
    main()
