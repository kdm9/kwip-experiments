#!/usr/bin/env python2

# Use either ete2 or ete3
try:
    import ete3 as ete
except ImportError:
    import ete2 as ete

import numpy as np

CLI = """
USAGE:
    random_subtree <tree> <n>

Subsamples <n> taxa from the Newick tree in <tree>, preserving the branch
lengths of subsampled taxa.
"""

def main(treefile, n):
    n = int(n)
    tree = ete.Tree(treefile)
    leaves = tree.get_leaf_names()
    subsample = [leaves[i] for i in np.random.choice(n, size=len(tree))]
    tree.prune(subsample, preserve_branch_length=True)
    print(tree.write())

if __name__ == "__main__":
    import docopt
    opts = docopt.docopt(CLI)
    main(opts['<tree>'], int(opts['<n>']))
