import numpy as np
import pandas as pd
import os
from tqdm import tqdm

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


# Function to generate random memberships
def synthetic_membership(genome_genes, target_genes):
    # Remove duplicates by converting lists to sets and back to lists
    genome_genes = list(set(genome_genes))
    target_genes = list(set(target_genes))
    
    if len(genome_genes) < 20:
        raise ValueError("Not enough genes to sample random genes. Provide a larger genome list.")
    
    # Ensure target_genes does not exceed genome_genes
    target_genes = [gene for gene in target_genes if gene in genome_genes]

    random_genes = np.random.choice(
        [gene for gene in genome_genes if gene not in target_genes], 
        size=int(0.1 * len(target_genes) / (1 - 0.1)), 
        replace=False
    )
    
    query_genes = set(target_genes).union(set(random_genes))
    query_membership = [1 if gene in query_genes else 0 for gene in genome_genes]
    
    return query_membership

def save_synthetic_memberships(n_it, genome_path, pathway_file, pathway_ids, output_dir):
    genome_genes = pd.read_csv(genome_path, sep='\t')['ENSEMBL'].tolist()
    genome_genes = list(set(genome_genes))  # Remove duplicates
    print(f"Unique genome genes: {len(genome_genes)}")  # Check unique genes

    os.makedirs(output_dir, exist_ok=True)

    for pathway_name in pathway_ids:
        target_genes = pd.read_csv(pathway_file, sep='\t')
        target_genes = target_genes[target_genes['Pathway_Name'] == pathway_name]['Ensembl_ID'].tolist()
        target_genes = list(set(target_genes))  # Remove duplicates
        print(f"Number of target genes for {pathway_name}: {len(target_genes)}")  # Check target genes

        random_query_memberships_list = []
        for _ in tqdm(range(n_it), desc=f"Generating synthetic memberships for {pathway_name}"):
            membership = synthetic_membership(genome_genes, target_genes)
            if len(membership) == len(genome_genes):
                random_query_memberships_list.append(membership)
            else:
                print(f"Generated membership of unexpected length {len(membership)} for pathway {pathway_name}.")

        # Construct DataFrame only if list lengths are consistent
        if random_query_memberships_list:
            random_query_memberships = pd.DataFrame(
                data=np.array(random_query_memberships_list).T,
                columns=[f"Iter_{i}_Membership" for i in range(n_it)]
            )
            
            # Check if the lengths match
            if len(random_query_memberships) == len(genome_genes):
                random_query_memberships['Ensembl_ID'] = genome_genes
            else:
                print(f"Length mismatch: {len(random_query_memberships)} vs {len(genome_genes)}")

            synthetic_membership_path = os.path.join(output_dir, f'synthetic_memberships_{pathway_name}.tsv')
            random_query_memberships.to_csv(synthetic_membership_path, sep='\t', index=False)
            print(f"Synthetic memberships for {pathway_name} saved to {synthetic_membership_path}")
        else:
            print(f"No valid memberships generated for pathway {pathway_name}.")

# Parameters
genome_path = '../../../data/id_conversion/ENTREZID_to_ENSEMBL.txt'
n_it = 10
pathway_file = "../../../data/pathways/KEGG/KEGG_2022_pathway_memberships.tsv"
output_dir = '../../../data/query/synthetic/'

# List of pathway IDs for which to generate synthetic memberships
pathway_ids = infection_pathway_ids 

# Run membership generation for each pathway ID
save_synthetic_memberships(n_it, genome_path, pathway_file, pathway_ids, output_dir)
