import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from joblib import Parallel, delayed, cpu_count

from fuzzy_gsea_t import fuzzy_gsea_score 


def t_stat_proc(data, labels):
    permuted_labels = np.random.permutation(labels)
    set1 = data.loc[:,labels]
    set2 = data.loc[:,~labels]
    t_scores, _ = ttest_ind(set1,set2, axis=1, equal_var=False, nan_policy='omit')
    return t_scores

def fast_sample_label_permutation(data, labels, n_permutations = 1000, n_jobs = cpu_count()-1):
    labels = np.array(labels, dtype=bool)
    result = Parallel(n_jobs=n_jobs)(delayed(t_stat_proc)(data, labels) for i in range(n_permutations))
    result = pd.DataFrame(np.array(result).transpose(), index=data.index)
    return result

def gsea_proc(tests_statistic, pathways_dict):
    test_statistic = tests_statistic.sort_values(ascending = False)
    result, _ = fuzzy_gsea_score(ranked_list, pathways_dict)
    return result

def gsea_parallel(permutation_tests, pathways_dict):
    result = Parallel(n_jobs=n_jobs)(delayed(gsea_proc)(permutation_tests.iloc[:,i], labels) for i in range(permutation_tests.shape[1]))
    return result
