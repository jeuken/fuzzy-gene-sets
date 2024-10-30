import numpy as np
import pandas as pd
import os
import argparse
import tqdm
from scipy.stats import hypergeom  
from typing import List, Optional


def load_query(query_file: str, query_membership_type: str) -> pd.DataFrame:
    '''Load and clean query data from a file.'''
    if not os.path.isfile(query_file):
        raise FileNotFoundError(f"Query file '{query_file}' does not exist.")
    
    query_df = pd.read_csv(query_file, sep='\t')
    
    if 'Ensembl_ID' not in query_df.columns or query_membership_type not in query_df.columns:
        raise ValueError(f"Query file must contain 'Ensembl_ID' and '{query_membership_type}' columns.")
    
    query_df = query_df.dropna().rename(columns={query_membership_type: 'Query_Membership'})
    
    return query_df

def load_pathways(pathway_file: str, pathway_membership_type: str) -> pd.DataFrame:
    '''Load and aggregate pathway data from a file.'''
    if not os.path.isfile(pathway_file):
        raise FileNotFoundError(f"Pathway file '{pathway_file}' does not exist.")
    
    pathway_df = pd.read_csv(
        pathway_file,
        sep='\t',
        usecols=['Pathway_Name', 'Description', 'Ensembl_ID', pathway_membership_type],
        dtype={'Ensembl_ID': str, pathway_membership_type: float}
    )
    
    if 'Pathway_Name' not in pathway_df.columns or 'Description' not in pathway_df.columns or pathway_membership_type not in pathway_df.columns:
        raise ValueError(f"Pathway file must contain 'Pathway_Name', 'Description', and '{pathway_membership_type}' columns.")
    
    pathway_df = (pathway_df
                   .dropna()
                   .rename(columns={pathway_membership_type: 'Pathway_Membership'})
                   .groupby('Pathway_Name')
                   .agg({
                       'Description': 'first',
                       'Ensembl_ID': list,
                       'Pathway_Membership': list
                   })
                   .reset_index())
    
    return pathway_df

def load_probabilities(probability_file: str) -> pd.DataFrame:
    '''Load precomputed probabilities from a file.'''
    if not os.path.isfile(probability_file):
        raise FileNotFoundError(f"Probability file '{probability_file}' does not exist.")
    
    probabilities_df = pd.read_csv(probability_file, sep='\t', index_col=0)
    return probabilities_df

def ora_fuzzy_intersection(query_memberships, pathway_memberships):
    '''Calculate the fuzzy intersection size between query and pathway memberships.'''
    intersection = np.multiply(query_memberships, pathway_memberships)
    intersection_size = np.sum(intersection)
    return intersection_size

def dp_p_value(observed_intersection, probabilities_df, universe_size, pathway_size, query_size, factor):
    '''Calculate the p-value based on the observed intersection using precomputed probabilities.'''
    # Divide the observed intersection by the user-defined factor
    adjusted_intersection = observed_intersection * factor

    # Round the adjusted intersection to 10 decimal places
    adjusted_intersection = round(adjusted_intersection, 10)

    try:
        closest_score = probabilities_df.index[probabilities_df.index.get_loc(adjusted_intersection)]  # Get closest score
    except KeyError:
        # If the adjusted intersection is not found, find the nearest index
        closest_score = probabilities_df.index[(np.abs(probabilities_df.index - adjusted_intersection)).argmin()]

    score_row = probabilities_df.loc[closest_score, :]  # Retrieve the row for that score

    if score_row.empty:
        raise ValueError(f"No entry found for adjusted intersection size: {adjusted_intersection}")

    p_value = 0
    for index, value in score_row.items():  
        try:
            index = int(index)  # Convert index to integer
        except ValueError:
            raise ValueError(f"Index {index} could not be converted to an integer.")

        prob = value * hypergeom.pmf(index, universe_size, pathway_size, query_size)
        p_value += prob

    return p_value 

def dp_ora_compute_stats(pathway, query_df, probability_file, factor, plots=False):
    '''Compute statistics for fuzzy ORA for a given pathway and query DataFrame.'''
    # Load pathway-specific probabilities file
    probabilities_df = load_probabilities(probability_file)

    pathway_df = pd.DataFrame({
        'Ensembl_ID': pathway['Ensembl_ID'],
        'Pathway_Membership': pathway['Pathway_Membership']
    })

    # Merge query_df with pathway_df
    merged_df = pd.merge(query_df, pathway_df, on='Ensembl_ID', how='left').fillna(0)
    
    query_memberships = merged_df['Query_Membership'].values
    pathway_memberships = merged_df['Pathway_Membership'].values
    
    observed_intersection = ora_fuzzy_intersection(query_memberships, pathway_memberships)

    # Calculate pathway size considering only pathway genes found in query_df
    pathway_size = merged_df[merged_df['Pathway_Membership'] > 0].shape[0]

    # Universe size is still based on merged_df
    universe_size = len(merged_df)
    # Query size is based on the number of genes with Query_Membership == 1
    query_size = len(merged_df[merged_df['Query_Membership'] == 1])
    
    p_value = dp_p_value(observed_intersection, probabilities_df, universe_size, pathway_size, query_size, factor)
    
    return observed_intersection, p_value

def dp_ora(
    query_file: str,
    pathway_file: str,
    probability_folder: str,
    factor: float = 113,  # User-defined factor with a default value
    query_membership_type: str = 'Crisp_Membership',
    pathway_membership_type: str = 'Overlap_Membership',
    output_path: Optional[str] = None,
    dataset_name: str = '',
    pathway_ids: Optional[List[str]] = None
) -> pd.DataFrame:
    '''Perform fuzzy Over-Representation Analysis (ORA) using precomputed probabilities for each pathway.'''
    
    query_df = load_query(query_file, query_membership_type)
    pathway_df = load_pathways(pathway_file, pathway_membership_type)

    if pathway_ids:
        pathway_df = pathway_df[pathway_df['Pathway_Name'].isin(pathway_ids)]

    results = []
    
    for idx, pathway in tqdm.tqdm(pathway_df.iterrows(), total=len(pathway_df), desc="Processing Pathways"):
        pathway_name = pathway['Pathway_Name']

        # Construct the file path for the pathway's specific probability file
        probability_file = os.path.join(probability_folder, f"{pathway_name}_dp_probabilities.tsv")
        
        # Check if the probability file exists
        if not os.path.exists(probability_file):
            print(f"Warning: Probability file for pathway '{pathway_name}' not found: {probability_file}. Skipping this pathway.")
            continue  # Skip to the next pathway

        observed_intersection, p_value = dp_ora_compute_stats(
            pathway, query_df, probability_file, factor  # Pass the factor here
        )

        results.append({
            'Pathway_Name': pathway_name,
            'Description': pathway['Description'],
            'Observed_Intersection': observed_intersection,
            'p_value': p_value
        })

    results_df = pd.DataFrame(results)
    results_df['p_value'] = results_df['p_value'].apply(lambda x: f"{x:.4e}")

    if output_path:
        output_file = os.path.join(output_path, f'{dataset_name}_{query_membership_type}_{pathway_membership_type}_results.csv')
        results_df.to_csv(output_file, index=False)
    
    return results_df


def dp_ora_main():
    """Parse command-line arguments and execute the dp_ora function."""
    parser = argparse.ArgumentParser(description="Run fuzzy ORA analysis.")
    parser.add_argument('-q', '--query_file', required=True, help="Path to the query file.")
    parser.add_argument('-p', '--pathway_file', required=True, help="Path to the pathway file.")
    parser.add_argument('-prob_folder', '--probability_folder', required=True, help="Path to the folder containing pathway-specific probability files.")
    parser.add_argument('-f', '--factor', type=float, default=113, help="Factor to adjust the observed intersection.")
    parser.add_argument('-q_name', '--query_membership_type', default='Crisp_Membership', help="Query membership type.")
    parser.add_argument('-p_name', '--pathway_membership_type', default='Crisp_Membership', help="Pathway membership type.")
    parser.add_argument('-o', '--output_path', default=None, help="Output directory path.")
    parser.add_argument('-d', '--dataset_name', default='dataset', help="Dataset name for output files.")
    parser.add_argument('-p_ids', '--pathway_ids', nargs='*', help="Specific pathway IDs to analyze.")

    args = parser.parse_args()
    
    dp_ora(
        query_file=args.query_file,
        pathway_file=args.pathway_file,
        probability_folder=args.probability_folder,
        factor=args.factor,  # Pass the user-defined factor
        query_membership_type=args.query_membership_type,
        pathway_membership_type=args.pathway_membership_type,
        output_path=args.output_path,
        dataset_name=args.dataset_name,
        pathway_ids=args.pathway_ids
    )


if __name__ == "__main__":
    # Uncomment the following line to run from the command line
    # dp_ora_main()

    # Direct invocation for IDE usage
    query_file = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data_old/single_cell/HIV_E-GEOD-111727_membership.csv"
    pathway_file =  '../../../data/pathways/KEGG/KEGG_2022_pathway_memberships.tsv'
    probability_folder = "../../../data/pathways/KEGG/DP/Overlap"
    query_membership_type = 'Crisp_Membership'
    pathway_membership_type = 'Overlap_Membership'
    output_path = "../../../data/ora_output/dp_ora_output"
    dataset_name = "Testing_DP"
    pathway_ids = None
    factor = 337

    dp_ora_results = dp_ora(
        query_file=query_file,
        pathway_file=pathway_file,
        probability_folder=probability_folder,
        query_membership_type=query_membership_type,
        pathway_membership_type=pathway_membership_type,
        output_path=output_path,
        dataset_name=dataset_name,
        pathway_ids=pathway_ids,
        factor = factor
    )

    print(dp_ora_results)

