from joblib import Parallel, delayed
from tqdm import tqdm
import os
import numpy as np
import pandas as pd
from scipy.stats import hypergeom, rankdata
import itertools


def load_pathway_data(pathway_file, membership_types):
    """Load pathway data and return a dictionary with pathway-specific DataFrames."""
    if not os.path.exists(pathway_file):
        raise FileNotFoundError(f"Pathway file '{pathway_file}' does not exist.")
    
    pathway_data = pd.read_csv(pathway_file, sep='\t', dtype={'Ensembl_ID': str})
    relevant_columns = ['Pathway_Name', 'Description', 'Ensembl_ID'] + membership_types
    pathway_data = pathway_data[relevant_columns].dropna(subset=['Pathway_Name', 'Ensembl_ID'])
    
    return {name: group.reset_index(drop=True) for name, group in pathway_data.groupby('Pathway_Name')}


def load_query(query_file, query_membership_type):
    """Load and clean query data."""
    if not os.path.isfile(query_file):
        raise FileNotFoundError(f"Query file '{query_file}' does not exist.")
    
    query_df = pd.read_csv(query_file, sep='\t')
    if 'Ensembl_ID' not in query_df.columns or query_membership_type not in query_df.columns:
        raise ValueError(f"Query file must contain 'Ensembl_ID' and '{query_membership_type}' columns.")
    
    return query_df.rename(columns={query_membership_type: 'Query_Membership'}).dropna()


def ora_fuzzy_intersection(query_memberships, pathway_memberships):
    """Calculate the fuzzy intersection size between query and pathway memberships."""
    return np.sum(np.multiply(query_memberships, pathway_memberships))


def calculate_hypergeometric_params(pathway_df, query_df, membership_type):
    """Calculate hypergeometric parameters for a given pathway and query."""
    
    # Prepare pathway dataframe with Ensembl_ID and Pathway Membership
    pathway_df = pd.DataFrame({
        'Ensembl_ID': pathway_df['Ensembl_ID'],
        'Pathway_Membership': pathway_df[membership_type]
    })
    
    # Merge query dataframe with pathway dataframe on Ensembl_ID
    merged_df = pd.merge(query_df, pathway_df, on='Ensembl_ID', how='left').fillna(0)
    
    # Extract query and pathway membership values
    query_memberships = merged_df['Query_Membership'].values
    pathway_memberships = merged_df['Pathway_Membership'].values
    
    # Calculate observed intersection: genes with non-zero pathway membership and crisp query membership of 1
    crisp_intersection = ((pathway_memberships > 0) & (query_memberships == 1)).sum()
    fuzzy_intersection = ora_fuzzy_intersection(query_memberships, pathway_memberships)
    
    # Calculate hypergeometric parameters
    pathway_size = (merged_df['Pathway_Membership'] > 0).sum()
    universe_size = len(merged_df)
    query_size = (merged_df['Query_Membership'] == 1).sum()
    
    return universe_size, query_size, pathway_size, crisp_intersection, fuzzy_intersection


def compute_score_distribution(prob_df, universe_size, query_size, pathway_size):
    # Step 1: Identify numeric columns in the DataFrame.
    numeric_columns = prob_df.columns[prob_df.columns.str.isnumeric()]
    if numeric_columns.empty:
        raise ValueError("No numeric columns found in the input DataFrame.")
    
    # Step 2: Convert numeric column names to floats to use as observed intersection values (k).
    columns = numeric_columns.astype(float)
    
    # Step 3: Calculate hypergeometric PMF for each possible score (k values).
    pmf_values = hypergeom.pmf(columns.astype(int), M=universe_size, n=pathway_size, N=query_size)
    
    # Check if the pmf_values sum to 1, else raise an error
    if not np.isclose(pmf_values.sum(), 1):
        raise ValueError("The hypergeometric PMF values do not sum to 1.")
    
    if pmf_values.size == 0:
        raise ValueError("Calculated PMF values are empty. Check input parameters.")
    
    # Step 4: Compute weighted probabilities by multiplying DataFrame values by PMF values.
    probabilities = (prob_df[numeric_columns].values * pmf_values).sum(axis=1)
    # Check if the probabilities sum to 1, else raise an error
    if not np.isclose(probabilities.sum(), 1):
        raise ValueError("The calculated probabilities do not sum to 1.")
    
    # Step 5: Return a DataFrame with the scores and their probabilities.
    score_dist = probabilities
    
    return score_dist

def create_probability_array(score_dist, factor):
    """Transform score distributions into an array by scaling indices."""
    if factor == 0:
        # If factor is zero, return None (no probability array for this membership)
        return None 
    
    # Initialize the maximum index in the new vector based on the original vector length
    max_index = len(score_dist) * factor * 10

    # Create a probability vector with zeros
    prob_vector = np.zeros(int(max_index) + 1)
    
    # Map each probability from score_dist to the scaled index in prob_vector
    for idx, prob in enumerate(score_dist):
        new_index = int(idx * factor * 10)
        prob_vector[new_index] = prob
    
    return prob_vector


def convolve_probabilities(prob_vectors):
    """Convolve multiple probability arrays."""
    combined_probs = prob_vectors[0]
    for prob_vector in prob_vectors[1:]:
        combined_probs = np.convolve(combined_probs, prob_vector)
    return combined_probs


def calculate_p_value(observed_score, combined_probs):
    """Calculate the p-value for an observed score."""
    return np.sum(combined_probs[observed_score:])


def process_pathway(pathway_name, pathway_df, query_df, probability_dir, membership_types):
    """Process a single pathway and calculate p-values for each parameter combination."""
    score_dists = {}
    observed_scores = {}

    # Load score distributions for each membership type
    for membership_type in membership_types:
        prob_file = os.path.join(probability_dir, membership_type, f"{pathway_name}_dp_probabilities.tsv")
        if os.path.exists(prob_file):
            prob_df = pd.read_csv(prob_file, sep="\t")
            universe_size, query_size, pathway_size, crisp_intersection, observed_score = calculate_hypergeometric_params(
                pathway_df, query_df, membership_type)
            observed_scores[membership_type] = observed_score
            score_dists[membership_type] = compute_score_distribution(
                prob_df, universe_size, query_size, pathway_size)
        else:
            return pathway_name, None

    # Parameter combinations for the fuzzy GSEA
    param_combinations = [
        params for params in itertools.product([0.0, 0.2, 0.4, 0.6, 0.8,1.0], repeat=len(membership_types))
        if sum(params) == 1
    ]
    pathway_p_values = []

    for params in param_combinations:
        prob_factors = []
        
        # Only add non-empty probability arrays (skip factor == 0)
        for m, param in zip(membership_types, params):
            prob_array = create_probability_array(score_dists[m], param)
            if prob_array is not None:  # Skip if the probability array is None (factor == 0)
                prob_factors.append(prob_array)
        
        # Skip the combination if no valid probability arrays (i.e., all factors are zero)
        if not prob_factors:
            continue
        
        # If there is only one valid probability vector, directly use it to calculate the p-value
        if len(prob_factors) == 1:
            combined_probs = prob_factors[0]
            score = sum(observed_scores[m] * param * 100 for m, param in zip(membership_types, params))
            p_val = calculate_p_value(int(score), combined_probs)
        else:
            # Perform convolution if multiple probability vectors are present
            combined_probs = convolve_probabilities(prob_factors)
            score = sum(observed_scores[m] * param * 100 for m, param in zip(membership_types, params))
            p_val = calculate_p_value(int(score), combined_probs)
        
        pathway_p_values.append((params, p_val))

    return pathway_name, pathway_p_values

def process_dataset(dataset, pathway_data, probability_dir, membership_types):
    """Process a single dataset and return the computed results."""
    condition = dataset["Disease"]
    target_pathway = dataset["Target Pathway"]
    dataset_name = dataset["Dataset"]
    query_file = f"../../../data/query/real/bulk/{condition}_{dataset_name}_membership_1_fc.csv"
    query_df = load_query(query_file, "Crisp_Membership")

    pathway_p_values = {}

    results = Parallel(n_jobs=-1)(
        delayed(process_pathway)(
            pathway_name, pathway_df, query_df, probability_dir, membership_types
        )
        for pathway_name, pathway_df in tqdm(pathway_data.items(), desc="Pathways", leave=False)
    )

    for pathway_name, p_values in results:
        if p_values:
            pathway_p_values[pathway_name] = p_values

    final_results = []
    for params, _ in pathway_p_values.get(target_pathway, []):
        sorted_pathways = sorted(
            [(name, p_val) for name, values in pathway_p_values.items()
             for param, p_val in values if param == params],
            key=lambda x: x[1]
        )
        if target_pathway not in [name for name, _ in sorted_pathways]:
            continue
        ranks = rankdata([p_val for _, p_val in sorted_pathways], method="average")
        target_rank = ranks[[name for name, _ in sorted_pathways].index(target_pathway)]
        
        # Get the p-value for the target pathway from the pathway_p_values
        target_p_value = next((p_val for param, p_val in pathway_p_values.get(target_pathway, [])
                               if param == params), None)
        
        # Create a result dictionary with separate columns for each membership parameter
        result = {
            "Dataset": dataset_name,
            "Target_Rank": int(target_rank),
            "Target_Pathway": target_pathway,
            "Disease": condition,
            "Target_P_Value": target_p_value
        }
        
        # Dynamically add parameter columns for each membership type
        for membership_type, param_value in zip(membership_types, params):
            result[f"{membership_type}_Parameter"] = param_value

        final_results.append(result)
    return final_results


def main(pathway_file, datasets_file, probability_dir, membership_types, output_dir,name):
    """Main function to process datasets and save results."""
    pathway_data = load_pathway_data(pathway_file, membership_types)
    datasets = pd.read_csv(datasets_file)

    results = Parallel(n_jobs=-1)(
        delayed(process_dataset)(
            dataset, pathway_data, probability_dir, membership_types
        )
        for dataset in tqdm(datasets.to_dict(orient="records"), desc="Datasets")
    )

    flattened_results = [result for sublist in results for result in sublist]

    output_file = os.path.join(output_dir, f"parameter_search_results_{name}.csv")
    pd.DataFrame(flattened_results).to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")


if __name__ == "__main__":
    pathway_file = "../../../data/pathways/KEGG/KEGG_2022_pathway_memberships.tsv"
    datasets_file = "../../../data/query/real/bulk/datasets/cancer_infection_datasets.csv"
    membership_types = ["Crisp_Membership", "Strict_Overlap_Membership", "Harmonic_In_Membership","Harmonic_Out_Membership"]
    probability_dir = "../../../data/pathways/KEGG/DP_strict/"
    output_dir = "../../../data/ora_output/parameter_search/"
    name = "cancer_infection_fc_crisp_strict_overlap_harmonic_in_harmonic_out"
    main(pathway_file, datasets_file, probability_dir, membership_types, output_dir,name)
