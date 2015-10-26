#!/usr/bin/env python
from __future__ import print_function
from skbio.stats.distance import DissimilarityMatrix
from os import path
import sys

for distfile in sys.argv[1:]:
    prefix = path.splitext(distfile)[0]
    print(prefix)
    dist = DissimilarityMatrix.read(distfile)
    fig = dist.plot(title=prefix)
    fig.tight_layout()
    fig.savefig("{prefix}.pdf".format(prefix=prefix))
