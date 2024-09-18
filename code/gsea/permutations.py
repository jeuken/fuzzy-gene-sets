import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from joblib import Parallel, delayed, cpu_count

def t_stat_proc(data, labels):
    permuted_labels = np.random.permutation(labels)
    set1 = data.loc[:,labels]
    set2 = data.loc[:,~labels]
    t_scores, _ = ttest_ind(set1,set2, axis=1, equal_var=False, nan_policy='omit')
    return t_scores


def fast_sample_label_permutation(data, labels, n_permutations = 1000, n_jobs = cpu_count()-1):
    labels = np.array(labels, dtype=bool)
    result = Parallel(n_jobs=n_jobs)(delayed(t_stat_proc)(data, labels) for i in range(n_permutations))
    return result
