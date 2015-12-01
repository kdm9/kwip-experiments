#!/usr/bin/env python
from __future__ import division, print_function
from cogent import LoadTree
import random
import string
from math import log, ceil
import itertools


CLI = """
USAGE:
    random_subtree [options] <tree>

OPTIONS:
    -n N    Number of taxa to sample.

Subsamples N taxa from the Newick tree in <tree>, preserving the branch
lengths of subsampled taxa.
"""


def iter_trees(treefile):
    with open(treefile) as tfh:
        for line in tfh:
            yield LoadTree(treestring=line)

def scaled_subsample(tree, num, seed=None):
    if seed:
        random.seed(seed)
    sample = random.sample(sorted(tree.getTipNames()), int(num))
    subtree = tree.getSubTree(sample)
    subtree.scaleBranchLengths(1000000, ultrametric=True)
    lablen = int(ceil(log(num, 26)))
    labels = [
        "".join(x)
        for x in itertools.product(string.ascii_uppercase, repeat=lablen)
    ]

    for i, tip in enumerate(sorted(subtree.iterTips(), key=lambda x: x.Name)):
        print(i, labels[i], tip.Name)
    return subtree.getNewick(with_distances=True)


def dawg_tree(newicks):
    """List of newick trees as strings => a DAWG tree file"""
    treestr = ',\n\t'.join(newicks)
    return "Tree = {\n\t%s,\n}\n" % treestr


def main(tree, n):
    seed = random.random()
    trees = [scaled_subsample(t, n, seed) for t in iter_trees(tree)]
    with open("data/tree.dawg", 'w') as dfh:
        print(dawg_tree(trees), file=dfh)


if __name__ == "__main__":
    import docopt
    opts = docopt.docopt(CLI)
    main(opts['<tree>'], int(opts['-n']))
