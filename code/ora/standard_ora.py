cancer_pathway_ids = [
    "hsa05210",  # Colorectal cancer
    "hsa05212",  # Pancreatic cancer
    "hsa05225",  # Hepatocellular carcinoma
    "hsa05226",  # Gastric cancer
    "hsa05214",  # Glioma
    "hsa05216",  # Thyroid cancer
    "hsa05221",  # Acute myeloid leukemia
    "hsa05220",  # Chronic myeloid leukemia
    "hsa05217",  # Basal cell carcinoma
    "hsa05218",  # Melanoma
    "hsa05211",  # Renal cell carcinoma
    "hsa05219",  # Bladder cancer
    "hsa05215",  # Prostate cancer
    "hsa05213",  # Endometrial cancer
    "hsa05224",  # Breast cancer
    "hsa05222",  # Small cell lung cancer
    "hsa05223"   # Non-small cell lung cancer
]

viral_pathway_ids = [
    "hsa05166",  # Human T-cell leukemia virus 1 infection
    "hsa05170",  # Human immunodeficiency virus 1 infection
    "hsa05161",  # Hepatitis B
    "hsa05160",  # Hepatitis C
    "hsa05171",  # Coronavirus disease - COVID-19
    "hsa05164",  # Influenza A
    "hsa05162",  # Measles
    "hsa05168",  # Herpes simplex virus 1 infection
    "hsa05163",  # Human cytomegalovirus infection
    "hsa05167",  # Kaposi sarcoma-associated herpesvirus infection
    "hsa05169",  # Epstein-Barr virus infection
    "hsa05165"   # Human papillomavirus infection
]

infection_pathways = [
    'hsa05166',  # Human T-cell leukemia virus 1 infection
    'hsa05170',  # Human immunodeficiency virus 1 infection
    'hsa05161',  # Hepatitis B
    'hsa05160',  # Hepatitis C
    'hsa05171',  # Coronavirus disease - COVID-19
    'hsa05164',  # Influenza A
    'hsa05162',  # Measles
    'hsa05168',  # Herpes simplex virus 1 infection
    'hsa05163',  # Human cytomegalovirus infection
    'hsa05167',  # Kaposi sarcoma-associated herpesvirus infection
    'hsa05169',  # Epstein-Barr virus infection
    'hsa05165',  # Human papillomavirus infection
    'hsa05110',  # Vibrio cholerae infection
    'hsa05120',  # Epithelial cell signaling in Helicobacter pylori infection
    'hsa05130',  # Pathogenic Escherichia coli infection
    'hsa05132',  # Salmonella infection
    'hsa05131',  # Shigellosis
    'hsa05135',  # Yersinia infection
    'hsa05133',  # Pertussis
    'hsa05134',  # Legionellosis
    'hsa05150',  # Staphylococcus aureus infection
    'hsa05152',  # Tuberculosis
    'hsa05146',  # Amoebiasis
    'hsa05144',  # Malaria
    'hsa05145',  # Toxoplasmosis
    'hsa05140',  # Leishmaniasis
    'hsa05142',  # Chagas disease
    'hsa05143'   # African trypanosomiasis
]



#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import scipy.stats as stats
import os
import argparse
import tqdm

def ora_get_query(query_file, query_membership_type):
    if not os.path.isfile(query_file):
        raise FileNotFoundError(f"The query file {query_file} does not exist.")
    
    query_df = pd.read_csv(query_file, sep='\t', usecols=['Ensembl_ID', query_membership_type])
    query_df.dropna(subset=[query_membership_type], inplace=True)
    query_df.rename(columns={query_membership_type: 'Membership'}, inplace=True)
    return query_df

def ora_get_pathways(pathway_file, pathway_membership_type):
    df = pd.read_csv(
        pathway_file, 
        sep='\t', 
        usecols=['Pathway_Name', 'Description', 'Ensembl_ID', pathway_membership_type],
        dtype={'Ensembl_ID': str, pathway_membership_type: float}
    ).dropna()
    df.rename(columns={pathway_membership_type: 'Membership'}, inplace=True)
    pathway_df = df.groupby('Pathway_Name').agg({
        'Description': 'first',
        'Ensembl_ID': list,
        'Membership': list
    }).reset_index()
    return pathway_df

def ora_hypergeometric_test(query_genes, pathway_genes, universe_size):
    overlap_size = len(set(query_genes).intersection(set(pathway_genes)))
    M = universe_size  # Total number of genes in the universe
    n = len(query_genes)  # Total number of query genes with Membership 1
    N = len(pathway_genes)  # Total number of genes in the pathway
    k = overlap_size  # Number of genes in both the query and the pathway (overlap)

    # Hypergeometric test 
    p_value = stats.hypergeom.sf(k - 1, M, N, n)
    
    return overlap_size, p_value

def standard_ora_compute_stats(pathway, query_df, universe_size):
    # Filter query genes based on Membership and create a list
    query_genes = query_df['Ensembl_ID'].tolist()
    query_genes_filtered = query_df[query_df['Membership'] == 1]['Ensembl_ID'].tolist()
    
    # Access the Ensembl_IDs and Memberships directly
    pathway_genes = pathway['Ensembl_ID']  # List of Ensembl IDs
    pathway_memberships = pathway['Membership']  # List of Memberships

    # Filter pathway genes based on Membership and that are in query_genes
    pathway_genes_filtered = [
        gene for gene, membership in zip(pathway_genes, pathway_memberships) 
        if membership == 1 and gene in query_genes  # Membership must be 1 and in filtered query_genes
    ]

    overlap_size, p_value = ora_hypergeometric_test(query_genes_filtered, pathway_genes_filtered, universe_size)
    
    return overlap_size, p_value

def standard_ora(query_file, pathway_file, query_membership_type='Crisp', pathway_membership_type='Crisp', output_path=None, dataset_name='dataset', pathway_ids=None):
    query_df = ora_get_query(query_file, query_membership_type)
    pathway_df = ora_get_pathways(pathway_file, pathway_membership_type)
    
    # Filter pathways if a list of pathway IDs is provided
    if pathway_ids:
        pathway_df = pathway_df[pathway_df['Pathway_Name'].isin(pathway_ids)]
    
    universe_size = len(query_df)  # Universe size is the total number of genes in the query file

    results = []
    num_pathways = len(pathway_df)
    print(f"Total pathways: {num_pathways}")

    for idx, pathway in tqdm.tqdm(pathway_df.iterrows(), total=num_pathways, desc="Processing Pathways"):
        overlap_size, p_value = standard_ora_compute_stats(pathway, query_df, universe_size)
        
        results.append({
            'Pathway_Name': pathway['Pathway_Name'],
            'Description': pathway['Description'],
            'Observed_Intersection': overlap_size,
            'p_value': p_value
        })

    results_df = pd.DataFrame(results)
    results_df.sort_values('p_value', inplace=True)
    results_df.reset_index(drop=True, inplace=True)
    results_df['Rank'] = results_df['p_value'].rank(method='min').astype(int)
    
    if output_path:
        results_folder = os.path.join(output_path, dataset_name)
        os.makedirs(results_folder, exist_ok=True)  # Create output directory if it doesn't exist
        results_file_name = f"standard_ora_results_{dataset_name}_{query_membership_type}_{pathway_membership_type}.csv"
        results_df.to_csv(os.path.join(results_folder, results_file_name), index=False)
    
    return results_df

def main():
    parser = argparse.ArgumentParser(description="Run standard ORA analysis.")
    parser.add_argument('-q', '--query_file', required=True, help="Path to the query file.")
    parser.add_argument('-p', '--pathway_file', required=True, help="Path to the pathway file.")
    parser.add_argument('-q_name', '--query_membership_type', default='Crisp_Membership', help="Query membership type.")
    parser.add_argument('-p_name', '--pathway_membership_type', default='Crisp_Membership', help="Pathway membership type.")
    parser.add_argument('-o', '--output_path', default=None, help="Output directory path.")
    parser.add_argument('-d', '--dataset_name', default='dataset', help="Dataset name for output files.")
    parser.add_argument('--pathway_ids', nargs='+', help="List of specific pathway IDs to include.")

    args = parser.parse_args()
    
    results_df = standard_ora(
        query_file=args.query_file,
        pathway_file=args.pathway_file,
        query_membership_type=args.query_membership_type,
        pathway_membership_type=args.pathway_membership_type,
        output_path=args.output_path,
        dataset_name=args.dataset_name,
        pathway_ids=args.pathway_ids  # Pass the list of pathway IDs
    )
    
    print("Results:\n", results_df.head())

if __name__ == "__main__":
    # Uncomment the following lines to run the script directly from an IDE.
    query_file = "../../data/single_cell/HIV_E-GEOD-111727_membership.csv"
    pathway_file = "../../data/pathways/KEGG/KEGG_2022_membership.tsv"
    query_membership_type = 'Crisp_Membership'
    pathway_membership_type = 'Crisp_Membership'
    output_path = "../../data/ORA_output/infection/HIV"
    dataset_name = "HIV_E-GEOD-111727"
    pathway_ids = infection_pathways
    results_df = standard_ora(
        query_file=query_file,
        pathway_file=pathway_file,
        query_membership_type=query_membership_type,
        pathway_membership_type=pathway_membership_type,
        output_path=output_path,
        dataset_name=dataset_name,
        pathway_ids=pathway_ids  # Pass the list of selected pathway IDs
    )
    print("Results:\n", results_df.head())

