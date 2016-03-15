#!/usr/bin/env python3
import skbio
from skbio import TreeNode, DistanceMatrix
from scipy import stats
import numpy as np
import pandas as pd
from os import path



def get_truth(treefile, reps):
    samp = skbio.read(treefile, into=TreeNode)
    dist = samp.tip_tip_distances(sorted(x.name for x in samp.tips()))
    runs = DistanceMatrix(
        np.repeat(np.repeat(dist.data, reps, axis=1), reps, axis=0))
    runs.ids = ['{}-{}'.format(g, i+1) for g in dist.ids for i in range(reps)]
    truth = runs.condensed_form()
    return truth


def get_spearmans(distfile, truth):
    distmat = DistanceMatrix.read(distfile)
    dist  = distmat.condensed_form()
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



