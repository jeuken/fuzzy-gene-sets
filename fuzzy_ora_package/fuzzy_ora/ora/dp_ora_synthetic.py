import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm

def load_pathway_ids(file_path):
    ids = {'cancer': [], 'infection': []}
    with open(file_path, 'r') as file:
        current_category = None
        for line in file:
            line = line.strip()
            if line.startswith('#'):
                if 'cancer' in line:
                    current_category = 'cancer'
                elif 'infection' in line:
                    current_category = 'infection'
            elif line:
                ids[current_category].append(line)
    return ids

from ora.dp_ora import dp_ora

def dp_ora_null(random_query_memberships_path, pathway_file, n_it, null_distribution_path, pathway_ids, probability_folder, pathway_membership_type, factor):
    pathways = pd.read_csv(pathway_file, sep='\t')['Pathway_Name'].unique()
    results = pd.DataFrame(columns=[f"Iter_{i}_p_value" for i in range(n_it)])

    # Filter only for specified pathway_ids
    pathways_to_include = [pathway for pathway in pathways if pathway in pathway_ids]
    
    for pathway in pathways_to_include:
        results.loc[pathway] = [None] * n_it  # Initialize rows for the specified pathways

    def run_iteration(i):
        column_name = f"Iter_{i}_Membership"
        iteration_result = dp_ora(
            query_file=random_query_memberships_path,
            pathway_file=pathway_file,
            query_membership_type=column_name,
            pathway_membership_type=pathway_membership_type,
            output_path=None,
            dataset_name='dataset',
            pathway_ids=pathway_ids,
            probability_folder=probability_folder,
            factor=factor
        )
        return iteration_result[['Pathway_Name', 'p_value']]

    iteration_results = Parallel(n_jobs=-1)(
        delayed(run_iteration)(i) for i in tqdm(range(n_it), desc="Running ORA iterations")
    )

    for i, iteration_result in enumerate(iteration_results):
        column_name = f"Iter_{i}_p_value"
        for pathway_name in iteration_result['Pathway_Name']:
            if pathway_name in results.index:
                results.at[pathway_name, column_name] = iteration_result.loc[
                    iteration_result['Pathway_Name'] == pathway_name, 'p_value'
                ].values[0]

    results.to_csv(null_distribution_path, sep='\t')
    print(f"Null distributions saved to {null_distribution_path}")


# Parameters
if __name__ == "__main__":
    n_it = 10
    pathway_file = "../../../data/pathways/KEGG/KEGG_2022_pathway_memberships.tsv"
    
    pathway_ids = load_pathway_ids("../../../data/pathways/KEGG/pathway_groups.txt")['infection']
    membership_type = 'Strict_Overlap_Membership'
    
    # Define probability folder based on membership type
    probability_folder = f"../../../data/pathways/KEGG/DP/{membership_type}"

    # Specify the factor
    factor = 113  # Change this value as needed

    # Loop through each infection pathway ID
    for pathway_id in pathway_ids:
        # Update synthetic membership path for each pathway_id
        synthetic_membership_path = f'../../../data/query/synthetic/synthetic_memberships_{pathway_id}.tsv'
        null_distribution_path = f'../../../data/query/synthetic/synthetic_null_distributions_dp_{membership_type.lower()}_{pathway_id}.csv'
        
        dp_ora_null(synthetic_membership_path, pathway_file, n_it, null_distribution_path, pathway_ids, probability_folder, membership_type, factor)
