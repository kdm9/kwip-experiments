#!/usr/bin/env python3
from collections import defaultdict
import os
from os import path
import docopt
from sys import stdout, stderr


CLI = """
USAGE:
    mash2kwipdist MASHDIST
"""


def fname2id(fname):
    fname = path.basename(fname)
    exts = ["_il.fastq.gz", ".fastq.gz"]
    for ext in exts:
        if fname.endswith(ext):
            fname = fname[:-len(ext)]
            break
    return fname


def convert(infile, outfile='/dev/stdout'):
    dists = defaultdict(dict)
    with open(infile) as fh:
        for line in fh:
            dist = line.strip().split('\t')
            id1 = fname2id(dist[0])
            id2 = fname2id(dist[1])
            dist = float(dist[2])
            dists[id1][id2] = dist

    with open(outfile, 'w') as ofile:
        ids = [''] + list(sorted(dists.keys()))
        print(*ids, sep='\t', file=ofile)
        for id1, row in sorted(dists.items()):
            rowdists = [it[1] for it in sorted(row.items())]
            print(id1, *rowdists, sep='\t', file=ofile)


if __name__ == '__main__':
    opt = docopt.docopt(CLI)
    convert(opt['MASHDIST'])
