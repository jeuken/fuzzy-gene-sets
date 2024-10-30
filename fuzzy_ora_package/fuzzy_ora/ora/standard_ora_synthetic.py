import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm
import os
import sys

infection_pathway_ids = [
    'hsa05166', 'hsa05170', 'hsa05161', 'hsa05160', 'hsa05171',
    'hsa05164', 'hsa05162', 'hsa05168', 'hsa05163', 'hsa05167',
    'hsa05169', 'hsa05165', 'hsa05110', 'hsa05120', 'hsa05130',
    'hsa05132', 'hsa05131', 'hsa05135', 'hsa05133', 'hsa05134',
    'hsa05150', 'hsa05152', 'hsa05146', 'hsa05144', 'hsa05145',
    'hsa05140', 'hsa05142', 'hsa05143'
]


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

# Add the project root to sys.path
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(script_dir, '..'))
sys.path.insert(0, project_root)

from ora.standard_ora import standard_ora


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

# Run ORA on synthetic query sets to generate null distribution
def standard_ora_null(random_query_memberships_path, pathway_file, n_it, null_distribution_path, pathway_ids):
    """Run ORA on synthetic query sets and collect null distributions."""
    pathways = pd.read_csv(pathway_file, sep='\t')['Pathway_Name'].unique()
    results = pd.DataFrame(index=pathways, columns=[f"Iter_{i}_p_value" for i in range(n_it)])

    def run_iteration(i):
        """Perform standard ORA for one iteration."""
        column_name = f"Iter_{i}_Membership"
        iteration_result = standard_ora(
            query_file=random_query_memberships_path,
            pathway_file=pathway_file,
            query_membership_type=column_name,
            pathway_membership_type='Crisp_Membership',
            output_path=None,
            dataset_name='dataset',
            pathway_ids=pathway_ids
        )
        return iteration_result[['Pathway_Name', 'p_value']]

    # Run iterations in parallel and fill results
    iteration_results = Parallel(n_jobs=-1)(
        delayed(run_iteration)(i) for i in tqdm(range(n_it), desc="Running ORA iterations")
    )

    # Fill results DataFrame with p-values
    for i, iteration_result in enumerate(iteration_results):
        column_name = f"Iter_{i}_p_value"
        for pathway_name in iteration_result['Pathway_Name']:
            results.at[pathway_name, column_name] = iteration_result.loc[
                iteration_result['Pathway_Name'] == pathway_name, 'p_value'
            ].values[0]

    # Filter results to include only specified pathways
    results = results.loc[results.index.isin(pathway_ids)]

    # Save null distributions to file
    results.to_csv(null_distribution_path, sep='\t')
    print(f"Null distributions saved to {null_distribution_path}")
    return results

# Parameters
if __name__ == "__main__":
    n_it = 10
    pathway_file = "../../../data/pathways/KEGG/KEGG_2022_pathway_memberships.tsv"
    
    # Load only infection pathway IDs
    pathway_ids = load_pathway_ids("../../../data/pathways/KEGG/pathway_groups.txt")['cancer']

    # Loop through each infection pathway ID to generate null distributions
    for pathway_id in pathway_ids:
        # Define paths for synthetic memberships and null distribution files
        synthetic_membership_path = f'../../../data/query/synthetic/synthetic_memberships_{pathway_id}.tsv'
        null_distribution_path = f'../../../data/query/synthetic/synthetic_null_distributions_standard_{pathway_id}.csv'
        
        # Run ORA null distribution function for each pathway
        standard_ora_null(synthetic_membership_path, pathway_file, n_it, null_distribution_path, pathway_ids)
