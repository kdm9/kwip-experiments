#!/usr/bin/python3
from skbio import TreeNode, DistanceMatrix, TabularMSA, DNA
from scipy.spatial.distance import hamming

def aln_distmat(alignment, reps=3):
    '''Calculate pairwise distances from a MSA of genomes'''
    aln = TabularMSA.read(alignment, constructor=DNA)
    aln.reassign_index(minter="id")
    dist = DistanceMatrix.from_iterable([seq.values for seq in aln],
                                        metric=hamming, keys=aln.index)
    return dist

def alndist_namer(alnpath):
    for ext in ['.gz', '.fasta']:
        if alnpath.endswith(ext):
            alnpath = alnpath[:-len(ext)]
    alnpath += '.dist'
    return alnpath

if __name__ == '__main__':
    import docopt

    cli = '''
    USAGE:
        alndist.py ALNPATH
    '''
    opts = docopt.docopt(cli)

    aln = opts['ALNPATH']
    distmat = aln_distmat(aln)
    distmat.write(alndist_namer(aln))

