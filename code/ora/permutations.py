import numpy as np
from joblib import Parallel, delayed, cpu_count

def permutaion_proc(query_memberships, pathway_memberships):
    permuted_set = np.random.permutation(query_memberships)
    intersection = np.minimum(permuted_set,pathway_memberships)
    size = np.sum(intersection)
    return size


def fast_gene_label_permutation(query_memberships, pathway_memberships, n_permutations = 1000, n_jobs = cpu_count()-1):
    result = Parallel(n_jobs=n_jobs)(delayed(permutaion_proc)(query_memberships, pathway_memberships) for i in range(n_permutations))
    return result
