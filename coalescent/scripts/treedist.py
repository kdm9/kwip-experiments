#!/usr/bin/env python3
import skbio
from skbio import TreeNode, DistanceMatrix
from scipy import stats
import numpy as np
import pandas as pd

from os import path
import re

def iter_newick_partitoned(fname):
    with open(fname) as fh:
        for line in fh:
            m = re.match(r'\[(.*)\](\(.*;)', line)
            if m is None:
                # Assume it's just a normal newick tree
                yield 1, TreeNode.read([line])
            else:
                l, t = m.groups()
                yield int(float(l)), TreeNode.read([t])

def partition_weighted_distance(nwkfile):
    partdist = []
    partsum = 0
    totaldist = None
    tipnames = None
    for partlen, tree in iter_newick_partitoned(nwkfile):
        partsum += partlen
        if tipnames is None:
            try:
                tipnames = list(map(str, sorted(int(x.name) for x in tree.tips())))
            except ValueError:
                tipnames = list(sorted(x.name for x in tree.tips()))
        dist = tree.tip_tip_distances(tipnames).data
        if totaldist is None:
            totaldist = np.zeros_like(dist)
        partdist.append((partlen, dist))
    for partlen, dist in partdist:
        scale = partlen / partsum
        totaldist += dist * scale
    return DistanceMatrix(totaldist, ids=tipnames)


def get_truth(treefile, reps):
    dist = partition_weighted_distance(treefile)
    runs = DistanceMatrix(
        np.repeat(np.repeat(dist.data, reps, axis=1), reps, axis=0))
    runs.ids = ['{}-{}'.format(g, i+1) for g in dist.ids for i in range(reps)]
    return runs


def get_spearmans(distfile, truth):
    distmat = DistanceMatrix.read(distfile)
    ids = list(sorted(distmat.ids))
    distmat = distmat.filter(ids)
    dist  = distmat.condensed_form()
    truth = truth.condensed_form()
    sp = stats.spearmanr(truth, dist)
    return sp.correlation


def parse_stats(filename):
    filename = path.splitext(path.basename(filename))[0]
    cov, scale, measure = filename.split('-')
    cov = cov.rstrip('x')
    return {'coverage': float(cov),
            'scale': float(scale),
            'measure': measure}

def get_table(treefile, distfiles, reps=3):
    truth = get_truth(treefile, reps)
    runs = []
    for distfile in distfiles:
        run = parse_stats(distfile)
        run['spearman'] = get_spearmans(distfile, truth)
        runs.append(run)
    return pd.DataFrame(runs)


if __name__ == '__main__':
    import docopt

    cli = '''
    USAGE:
        treedist.py [-r REPS] -s SAMPLETREE <distmat> ...

    OPTIONS:
        -s SAMPLETREE       sample tree file.
        -r REPS             Reps per sample [default: 3]
    '''
    opts = docopt.docopt(cli)

    tree = opts['-s']
    reps = int(opts['-r'])
    distmats = opts['<distmat>']

    df = get_table(tree, distmats, reps)
    print(df.to_csv(index=False))
