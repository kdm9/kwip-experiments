#!/usr/bin/env python3
import numpy as np
from scipy import stats
from scipy.spatial.distance import hamming
from skbio import TreeNode, DistanceMatrix, TabularMSA, DNA
from docopt import docopt

import re


def sample_matrix_to_runs(dist, reps=3):
    '''Repeats a distance matrix to expand samples to reps.'''
    runs = DistanceMatrix(
        np.repeat(np.repeat(dist.data, reps, axis=1), reps, axis=0))
    runs.ids = ['{}-{}'.format(g, i+1) for g in dist.ids for i in range(reps)]
    return runs


def spearman(a, b):
    '''Returns spearman's \rho between a and b'''
    return stats.spearmanr(a, b).correlation


def distmat_corr(truthfile, distfile, reps=3, corrstat=spearman):
    '''Returns correlation between condensed distance matrices, using corrstat'''
    distmat = DistanceMatrix.read(distfile)
    truthmat = DistanceMatrix.read(truthfile)
    truthmat = sample_matrix_to_runs(truthmat, reps)

    ids = list(sorted(distmat.ids))
    t_ids = list(sorted(truthmat.ids))
    assert ids == t_ids, (ids, t_ids)

    dist = distmat.filter(ids).condensed_form()
    truth = truthmat.filter(ids).condensed_form()
    return corrstat(truth, dist)


CLI = '''
USAGE:
    calc_rho.py [options] TRUTH OBTAINED

OPTIONS:
  -S SEED    Seed
  -m METRIC  Metric
  -s SIZE    Sketchsize
  -c COV     Coverage
  -v VAR     Pi (mean pw dist)
  -r REPS    Replicate runs per sample
'''


def main():
    opts = docopt(CLI)

    rho = distmat_corr(opts["TRUTH"], opts["OBTAINED"], int(opts["-r"]))
    seed = opts['-S']
    metric = opts['-m']
    size = opts['-s']
    cov = opts['-c']
    var = opts['-v']
    print(seed, metric, size, cov, var, rho, sep='\t')

if __name__ == "__main__":
    main()
