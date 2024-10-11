import numpy as np
import pandas as pd
from joblib import Parallel, delayed, cpu_count
import matplotlib.pyplot as plt
import os
import argparse
import tqdm

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


top_infection_pathways = [
'hsa05169',
'hsa05168',
'hsa05166',
'hsa05164',
'hsa05170'
]



epstein = ['hsa05169']



def ora_get_query(query_file, query_membership_type):
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

def ora_fuzzy_intersection(query_memberships, pathway_memberships):
    intersection = np.multiply(query_memberships, pathway_memberships)  # Element-wise product
    intersection_size = np.sum(intersection)  # Sum of the product values
    return intersection_size


def ora_permutation(query_memberships, pathway_memberships):
    random_intersection_size = ora_fuzzy_intersection(np.random.permutation(query_memberships), pathway_memberships)
    return random_intersection_size

def ora_null_distribution(query_memberships, pathway_memberships, n_permutations=1000, n_jobs=cpu_count()-1):
    null_dist = Parallel(n_jobs=n_jobs)(delayed(ora_permutation)(query_memberships, pathway_memberships) for _ in range(n_permutations))
    return null_dist

def ora_p_value(observed_intersection, null_distribution):
    p_value = np.mean(null_distribution >= observed_intersection)
    return p_value

def ora_compute_stats(pathway, query_df, n_permutations, n_jobs, plots=False):
    pathway_df = pd.DataFrame({
        'Ensembl_ID': pathway['Ensembl_ID'],
        'Membership': pathway['Membership']
    })

    # Perform a left join to keep all query genes and include only matching pathway genes
    merged_df = pd.merge(query_df, pathway_df, on='Ensembl_ID', how='left', suffixes=('_query', '_pathway')).fillna(0)

    # Set pathway memberships to zero for genes not found in the pathway
    merged_df['Membership_pathway'] = merged_df['Membership_pathway'].fillna(0)

    # Get membership arrays for both query and pathway
    query_memberships_full = merged_df['Membership_query'].values
    pathway_memberships_full = merged_df['Membership_pathway'].values

    # Calculate observed intersection
    observed_intersection = ora_fuzzy_intersection(query_memberships_full, pathway_memberships_full)
    
    if plots:
        null_distribution = ora_null_distribution(query_memberships_full, pathway_memberships_full, n_permutations, n_jobs)
        p_value = ora_p_value(observed_intersection, null_distribution)
        return observed_intersection, null_distribution, p_value
    else:
        null_distribution = ora_null_distribution(query_memberships_full, pathway_memberships_full, n_permutations // 10, n_jobs)
        p_value = ora_p_value(observed_intersection, null_distribution)
        return observed_intersection, p_value

def ora_plot_null_distribution(pathway, observed_score, null_distribution, p_value, plot_path, query_membership_type, pathway_membership_type, dataset_name):
    plt.figure(figsize=(8, 6))
    plt.hist(null_distribution, bins=30, alpha=0.7, color='gray', edgecolor='black')
    plt.axvline(observed_score, color='red', linestyle='--', linewidth=2, label=f'Observed Score = {observed_score:.2f}')
    plt.title(f'Null Distribution for {pathway}\n(Query: {query_membership_type}, Pathway: {pathway_membership_type})')
    plt.xlabel('Fuzzy Intersection Size')
    plt.ylabel('Frequency')
    plt.annotate(f'P-value = {p_value:.4f}', xy=(0.7, 0.9), xycoords='axes fraction', fontsize=12, color='red')
    plt.legend()
    plt.tight_layout()
    
    folder_name = f"{query_membership_type}_{pathway_membership_type}"
    full_plot_path = os.path.join(plot_path, 'null_distributions_product', folder_name)
    os.makedirs(full_plot_path, exist_ok=True)
    
    file_name = f"{pathway}_{dataset_name}_{query_membership_type}_{pathway_membership_type}.png"
    plt.savefig(os.path.join(full_plot_path, file_name))
    plt.close()

def fuzzy_ora(query_file, pathway_file, query_membership_type='Crisp', pathway_membership_type='Crisp', n_permutations=1000, n_jobs=cpu_count()-1, output_path=None, plots=False, dataset_name='dataset', pathway_ids=None):
    query_df = ora_get_query(query_file, query_membership_type)
    pathway_df = ora_get_pathways(pathway_file, pathway_membership_type)
    
    if pathway_ids:
        # Filter pathway_df to include only pathways from the pathway_ids list
        pathway_df = pathway_df[pathway_df['Pathway_Name'].isin(pathway_ids)]
    
    results = []
    num_pathways = len(pathway_df)
    print(f"Total pathways: {num_pathways}")

    for idx, pathway in tqdm.tqdm(pathway_df.iterrows(), total=num_pathways, desc="Processing Pathways"):
        if plots:
            observed_intersection, null_distribution, p_value = ora_compute_stats(
                pathway, query_df, n_permutations, n_jobs, plots=True
            )
            plot_path = os.path.join(output_path, dataset_name) if output_path else os.path.join('.', dataset_name)
            ora_plot_null_distribution(
                pathway['Pathway_Name'], observed_intersection, null_distribution, p_value, 
                plot_path, query_membership_type, pathway_membership_type, dataset_name
            )
        else:
            observed_intersection, p_value = ora_compute_stats(
                pathway, query_df, n_permutations, n_jobs, plots=False
            )

        results.append({
            'Pathway_Name': pathway['Pathway_Name'],
            'Description': pathway['Description'],
            'Observed_Intersection': observed_intersection,
            'p_value': p_value
        })

    results_df = pd.DataFrame(results)
    results_df.sort_values('p_value', inplace=True)
    results_df.reset_index(drop=True, inplace=True)
    results_df['Rank'] = results_df['p_value'].rank(method='min').astype(int)
    
    if output_path:
        results_folder = os.path.join(output_path, dataset_name)
        os.makedirs(results_folder, exist_ok=True)
        results_file_name = f"ora_results_{dataset_name}_{query_membership_type}_{pathway_membership_type}.csv"
        results_df.to_csv(os.path.join(results_folder, results_file_name), index=False)
    
    return results_df



def main():
    parser = argparse.ArgumentParser(description="Run fuzzy ORA analysis.")
    parser.add_argument('-q', '--query_file', required=True, help="Path to the query file.")
    parser.add_argument('-p', '--pathway_file', required=True, help="Path to the pathway file.")
    parser.add_argument('-q_name', '--query_membership_type', default='Crisp_Membership', help="Query membership type.")
    parser.add_argument('-p_name', '--pathway_membership_type', default='Crisp_Membership', help="Pathway membership type.")
    parser.add_argument('-o', '--output_path', default=None, help="Output directory path.")
    parser.add_argument('-d', '--dataset_name', default='dataset', help="Dataset name for output files.")
    parser.add_argument('-n', '--n_permutations', type=int, default=1000, help="Number of permutations.")
    parser.add_argument('-j', '--n_jobs', type=int, default=cpu_count()-1, help="Number of parallel jobs.")
    parser.add_argument('--plots', action='store_true', help="Generate null distribution plots.")
    
    args = parser.parse_args()
    
    fuzzy_ora(
        query_file=args.query_file,
        pathway_file=args.pathway_file,
        query_membership_type=args.query_membership_type,
        pathway_membership_type=args.pathway_membership_type,
        n_permutations=args.n_permutations,
        n_jobs=args.n_jobs,
        output_path=args.output_path,
        plots=args.plots,
        dataset_name=args.dataset_name
    )

if __name__ == "__main__":
#    main()


# Comment the previous line and uncomment the following lines if you want to run this script directly from an IDE.
# This will bypass the command-line argument parsing and allow you to set parameters directly in the script.
    pathway_ids = viral_pathway_ids
    query_file = "../../data/single_cell/HIV_E-GEOD-111727_membership.csv"
    pathway_file = "../../data/pathways/KEGG/KEGG_2022_membership.tsv"
    query_membership_type = 'Crisp_Membership'
    pathway_membership_type = 'Overlap_Membership_linear'
    output_path = "../../data/ORA_output/infection/HIV"
    dataset_name = "HIV_E-GEOD-111727"
    pathway_ids = infection_pathways
    plots = False
    results_df = fuzzy_ora(
        query_file=query_file,
        pathway_file=pathway_file,
        query_membership_type=query_membership_type,
        pathway_membership_type=pathway_membership_type,
        n_permutations=10000000,
        n_jobs=cpu_count()-1,
        output_path=output_path,
        plots=False,
        dataset_name=dataset_name,
        pathway_ids=pathway_ids  # Pass the list of selected pathway IDs
    )

    print("Results:\n", results_df.head())



