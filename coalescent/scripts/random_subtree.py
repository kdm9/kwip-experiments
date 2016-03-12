#!/usr/bin/env python
from __future__ import division, print_function
from cogent import LoadTree
from skbio import TreeNode
import random
import string
from math import log, ceil
import itertools


CLI = """
USAGE:
    random_subtree [options] <tree>

OPTIONS:
    -o OUT      Output file.
    -n N        Number of taxa to sample.
    -s SEED     Seed [default: 1243]

Subsamples N taxa from the Newick tree in <tree>, preserving the branch
lengths of subsampled taxa.
"""


def iter_trees(treefile):
    with open(treefile) as tfh:
        for line in tfh:
            yield TreeNode
            (treestring=line)


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
        tip.Name = labels[i]
    return subtree.getNewick(with_distances=True)


def dawg_tree(newicks):
    """List of newick trees as strings => a DAWG tree file"""
    treestr = ',\n\t'.join(newicks)
    return "Tree = {\n\t%s,\n}\n" % treestr


def main(tree, n, out, seed):
    random.seed(seed)
    with open(out, 'w') as tfh:
        for tree in iter_trees(tree):
            subtree = scaled_subsample(t, n, seed=random.random()
            print(subtree, file=tfh)


if __name__ == "__main__":
    import docopt
    opts = docopt.docopt(CLI)
    main(opts['<tree>'], int(opts['-n']), opts['-o'], int(opts['-s']))
