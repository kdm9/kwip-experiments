import numpy as np
from scipy import stats
from scipy.spatial.distance import hamming
from skbio import TreeNode, DistanceMatrix, TabularMSA, DNA

import re


def aln_distmat(alignment, reps=3):
    '''Calculate pairwise distances from a MSA of genomes'''
    aln = TabularMSA.read(alignment, constructor=DNA)
    aln.reassign_index(minter="id")
    dist = DistanceMatrix.from_iterable([seq.values for seq in aln],
                                        metric=hamming, keys=aln.index)
    return dist


def iter_newick_partitoned(fname):
    '''Iterator over the trees in a partitioned newick file'''
    with open(fname) as fh:
        for line in fh:
            m = re.match(r'\[(.*)\](\(.*;)', line)
            if m is None:
                # Assume it's just a normal newick tree
                yield 1, TreeNode.read([line])
            else:
                l, t = m.groups()
                yield int(float(l)), TreeNode.read([t])


def tree_distmat(nwkfile):
    '''Calculate pairwise distance from a tree, weighted by any partitions'''
    partdist = []
    partsum = 0
    totaldist = None
    tipnames = None
    for partlen, tree in iter_newick_partitoned(nwkfile):
        partsum += partlen
        if tipnames is None:
            try:
                tipnames = list(map(
                    str, sorted(int(x.name) for x in tree.tips())))
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


def sample_matrix_to_runs(dist, reps=3):
    '''Repeats a distance matrix to expand samples to reps.'''
    runs = DistanceMatrix(
        np.repeat(np.repeat(dist.data, reps, axis=1), reps, axis=0))
    runs.ids = ['{}-{}'.format(g, i) for g in dist.ids for i in range(reps)]
    return runs


def load_sample_matrix_to_runs(samplematfile, reps=3):
    '''Loads a truth distance matrix between samples and expands to runs'''
    samples = DistanceMatrix.read(samplematfile)
    return sample_matrix_to_runs(samples, reps)


def spearmans_rho_distmats(distmat, truthmat):
    '''Calculates spearman's ρ between truth and dist's values. returns ρ'''
    ids = list(sorted(distmat.ids))
    t_ids = list(sorted(truthmat.ids))
    assert ids == t_ids, (ids, t_ids)
    dist = distmat.filter(ids).condensed_form()
    truth = truthmat.filter(ids).condensed_form()
    sp = stats.spearmanr(truth, dist)
    return sp.correlation
