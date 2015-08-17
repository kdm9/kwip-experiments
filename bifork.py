from __future__ import print_function
import gzip
from utils import *

genome = make_rand_genome(mbp=10)


bifork = biforcating_sequences(genome, levels=3, av_rate=0.0001, sd_rate=0.00001)
seqlist = list(flatten(bifork))

with open("data/bifork_10mb.fas", 'w') as fh:
    print_multifasta(seqlist, file=fh)


for i, seq in enumerate(seqlist):
    print('Genome', i)

    with gzip.open('data/gen_{}.fasta'.format(i), 'wb') as fh:
        fh.write(">{}\n{}\n".format(i, seq).encode('ascii'))

    n_reads = max(int(random.gauss(2e6, 4e5)), 10000)

