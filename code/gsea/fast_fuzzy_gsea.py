#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from joblib import Parallel, delayed, cpu_count


def gsea_get_pathways(pathway_file, pathway_membership_type):
    """Load pathway data and format it into a dictionary."""
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


def gsea_p_value(observed_score, null_distribution):
    """Calculate p-value for an observed score based on a null distribution."""
    null_distribution_abs = np.abs(null_distribution)
    observed_score_abs = abs(observed_score)
    p_value = np.mean(null_distribution_abs >= observed_score_abs)
    return p_value


def gsea_t_stat(data, labels):
    """Compute T-statistics for control and disease groups in expression data."""
    set_control = data.loc[:, labels == 0]  # control
    set_disease = data.loc[:, labels == 1]  # disease
    
    # Subtract control from disease to ensure positive t-stats when disease > control
    t_scores, _ = ttest_ind(set_disease, set_control, axis=1, equal_var=False, nan_policy='omit')
    return t_scores


def fast_sample_label_permutation(data, labels, n_permutations=1000, n_jobs=cpu_count() - 1):
    """Run permutations to generate a null distribution in parallel."""
    labels = np.array(labels, dtype=bool)

    def permute_labels():
        permuted_labels = np.random.permutation(labels)
        return gsea_t_stat(data, permuted_labels)

    results = Parallel(n_jobs=n_jobs)(delayed(permute_labels)() for _ in range(n_permutations))
    results_df = pd.DataFrame(np.array(results).T, index=data.index)
    return results_df


def gsea_proc(t_stats, pathways_dict):
    """Compute fuzzy GSEA scores for a set of T-statistics."""
    ranked_list = t_stats.sort_values(ascending=False).to_dict()
    enrichment_scores, _ = fuzzy_gsea_score(ranked_list, pathways_dict)
    return enrichment_scores


def gsea_parallel(permutation_tests, pathways_dict, n_jobs=cpu_count() - 1):
    """Parallelize GSEA score calculations for each permutation."""
    results = Parallel(n_jobs=n_jobs)(
        delayed(gsea_proc)(permutation_tests.iloc[:, i], pathways_dict) for i in range(permutation_tests.shape[1])
    )
    return results


def infer_labels_from_columns(expression_data):
    """Infer labels from column names (assuming 'control' and 'disease' are part of column names)."""
    column_names = expression_data.columns
    labels = np.where(column_names.str.contains('control', case=False, regex=True), 0, 1)  # 0 = control, 1 = disease
    return labels


def fuzzy_gsea(pathway_file, pathway_membership_type, expression_path, output_path, n_permutations=1000):
    """Main function to run fuzzy GSEA with parallelized permutations."""
    # Load pathways
    pathways = gsea_get_pathways(pathway_file, pathway_membership_type)
    
    # Load expression data
    expression_data = pd.read_csv(expression_path, index_col=0)

    # Infer labels from column names
    labels = infer_labels_from_columns(expression_data)

    # Compute original T-statistics
    original_t_stats = gsea_t_stat(expression_data, labels)

    # Generate null distribution through label permutation
    permuted_t_stats = fast_sample_label_permutation(expression_data, labels, n_permutations)

    # Compute enrichment scores for original data
    observed_scores = gsea_proc(original_t_stats, pathways)

    # Compute null distribution of enrichment scores in parallel
    null_distributions = gsea_parallel(permuted_t_stats, pathways)

    # Save observed scores and null distributions
    observed_scores.to_csv(f"{output_path}/observed_scores.csv", index=False)
    pd.DataFrame(null_distributions).to_csv(f"{output_path}/null_distributions.csv", index=False)

    print(f"Results saved to {output_path}")


if __name__ == "__main__":
    expression_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/gemma/Alzheimers_GSE95587_formatted.csv"
    pathway_file_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/pathways/KEGG/KEGG_2022_entrez_with_membership_new.tsv"
    output_path = "/Users/piamozdzanowski/VU/Fuzzy_Gene_Sets/data/GSEA_output/Alzheimers_GSE95587_trial"
    
    # Run the fuzzy GSEA analysis and save results
    fuzzy_gsea(pathway_file_path, 'Overlap_Membership_max', expression_path, output_path, n_permutations=1000)

