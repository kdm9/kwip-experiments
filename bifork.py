
# coding: utf-8

# In[ ]:

from utils import *


# Make a random genome, and 16 samples derived from it. Write it to a fasta file.

# In[ ]:

genome = make_rand_genome(mbp=10)


# In[ ]:

bifork = biforcating_sequences(genome, levels=3, av_rate=0.0001, sd_rate=0.00001)
seqlist = list(flatten(bifork))


# In[ ]:

with open("data/bifork_10mb.fas", 'w') as fh:
    print_multifasta(seqlist, file=fh)


# ## Make NJ tree
# 

# In[ ]:

from skbio import Alignment, DNA
from skbio.tree import nj


# In[ ]:

aln = Alignment.read('data/bifork_10mb.fas')


# In[ ]:

distmat = aln.distances()


# In[ ]:

distmat


# In[ ]:

tree = nj(distmat)


# In[ ]:

tree.write('data/bifork_10mb.nwk')


# In[ ]:

print(tree.ascii_art())


# ## Generate reads

# In[ ]:

from wgsim import *
import gzip
import random


# In[ ]:

for i, seq in enumerate(seqlist):
    print('Genome', i)
    
    with gzip.open('data/gen_{}.fasta'.format(i), 'wb') as fh:
        fh.write(">{}\n{}\n".format(i, seq).encode('ascii'))
    
    n_reads = max(int(random.gauss(2e6, 4e5)), 10000)
    get_ipython().system(' wgsim -1 101 -2 101 -r 0 -N {n_reads} data/gen_{i}.fasta /dev/shm/1.fq /dev/shm/2.fq')
    get_ipython().system(' pairs join  /dev/shm/1.fq /dev/shm/2.fq | gzip > "data/gen_{i}_il.fq.gz"')


# In[ ]:



