from utils import *
import sys
g = make_rand_genome(mbp=100)
print("made genome", file=sys.stderr)
b = biforcating_sequences(g, levels=4, av_rate=0.0001, sd_rate=0.00001)
print("made biforcations", file=sys.stderr)
sl = list(flatten(b))
print("flattened list", file=sys.stderr)
print(seqlist_to_mulitfasta(sl))
print("Done", file=sys.stderr)
