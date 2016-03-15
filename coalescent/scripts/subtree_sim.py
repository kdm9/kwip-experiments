#!/usr/bin/env python
import numpy as np
from cogent import LoadTree


CLI = """
USAGE:
    random_subtree <tree> <n>

Subsamples <n> taxa from the Newick tree in <tree>, preserving the branch
lengths of subsampled taxa.
"""

def main(treefile, n):
    n = int(n)
    tree = LoadTree(
    with open(treefile) as trees:
        for tree in trees:
            tree = ete.Tree(tree.strip())
            leaves = tree.get_leaf_names()
            subsample = [leaves[i] for i in np.random.choice(n, size=len(tree))]
            tree.prune(subsample, preserve_branch_length=True)
            print(tree.write())

if __name__ == "__main__":
    import docopt
    opts = docopt.docopt(CLI)
    main(opts['<tree>'], int(opts['<n>']))
