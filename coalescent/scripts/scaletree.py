#!/usr/bin/env python3
from __future__ import division, print_function

from ete3 import Tree
import numpy as np
from numpy import median, mean

import itertools as itl
from math import ceil, log
from string import ascii_uppercase


METRICS = [
    'mean',
    'median',
    'min',
    'max',
]


def pwdist(tree):
    '''Finds the (off-diagonal) pairwise distances between all tips of `tree`.
    '''
    dists = []
    for a, b in itl.combinations(tree.get_leaf_names(), 2):
        a = tree&a
        dists.append(a.get_distance(b))
    return np.array(dists)


def normalise_tree(tree, to=1.0, metric='mean'):
    '''
    Normalise branch lengths of `tree` such that `metric` of the pairwise
    distances is `to`.

    By default, normalise such that the mean of all pairwise distances is 1.0.
    '''
    dists = pwdist(tree)
    assert metric in METRICS
    current = eval('{}(dists)'.format(metric))

    for node in tree.iter_descendants():
        node.dist /= current
    return tree


def alphbetise_names(tree):
    '''Replace numeric tip labels with alphabetic ones. 1 -> A, 2 -> B etc.

    If there are more than 26 tips, labels are AA, AB, ..., ZZ and so forth for
    any number of tips.
    '''
    label_len = ceil(log(len(tree)) / log(26))  # how many letters do we need?
    labels = [''.join(letters)
              for letters in itl.product(ascii_uppercase, repeat=label_len)]

    tiplabels = list(sorted(tree.get_leaf_names(), key=int))
    for i, leaf in enumerate(tiplabels):
        node = tree&leaf
        node.name = labels[i]
    return tree


def main(treefile, to, metric):
    with open(treefile) as fh:
        for treeline in fh:
            tree = Tree(treeline)
            tree = alphbetise_names(tree)
            tree = normalise_tree(tree, to, metric)
            print(tree.write(format=5))


CLI = '''
USAGE:
    scaletree [options] TREEFILE

OPTIONS:
    -t TO       Scale tree metric to TO. [default: 1.0]
    -m METRIC   Metric for scaling. Must be one of mean, min, max.
                [default: mean]
'''

if __name__ == '__main__':
    from docopt import docopt
    opts = docopt(CLI)

    treefile = opts['TREEFILE']
    to = float(opts['-t'])
    metric = opts['-m']

    main(treefile, to, metric)
