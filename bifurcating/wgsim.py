import random
import sys
from math import log10
from utils import mutate_sequence

NT = list("ACGT")


def mutate_ref(ref, seed=None, mutrate=0.001, indelfrac=0.1, indelext=0.3,
               haploid=False):
    random.seed(seed)
    rand = random.random  # For brevity below
    haplotypes = [list(ref)]
    if not haploid:
        haplotypes.append(list(ref))
    idx = 0
    while idx < min(map(len, haplotypes)):
        if rand() < mutrate:
            homozygous = rand() < 0.333333
            strand = random.choice([0, 1])
            other = (strand + 1) % 2
            if rand() < indelfrac:
                deletion = rand() < 0.5
                while idx < min(map(len, haplotypes)):
                    if deletion:
                        haplotypes[strand].pop(idx)
                        if homozygous:
                            haplotypes[other].pop(idx)
                    else:
                        new_nt = random.choice(NT)
                        haplotypes[strand].insert(idx + 1, new_nt)
                        if homozygous:
                            haplotypes[other].insert(idx + 1, new_nt)
                    if rand() > indelext:
                        break
            else:
                old_nt = haplotypes[strand][idx]
                new_nt = old_nt
                while new_nt == old_nt:
                    new_nt = random.choice(NT)
                haplotypes[strand][idx] = new_nt
                if homozygous:
                    haplotypes[other][idx] = new_nt
        idx += 1
    return [''.join(h) for h in haplotypes]

C = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'a': 'T',
    'c': 'G',
    'g': 'C',
    't': 'A',
}

ILLUMINA = (
    'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT',
    'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT',
)

def revcomp(seq):
    new = []
    for x in reversed(seq):
        new.append(C.get(x, 'N'))
    return ''.join(new)

def fq(name, seq, qual):
    lines = [
        '@{}'.format(name),
        seq,
        '+',
        qual,
        '' # add a final newline
    ]
    return '\n'.join(lines)

def sim_reads(haplotypes, num_pairs, readlen=101, errrate=0.01, pair_dist=500,
              pair_sd=30, seed=None, adaptors=ILLUMINA, spos=None):
    random.seed(seed)

    qual_char = 'I'  # Maximum score in PHRED 33, for errrate == 0
    if errrate > 0:
        qual_char = chr(int(-10.0 * log10(errrate)) + 33)
    qual = qual_char * readlen

    for pair_num in range(num_pairs):
        this_pair_dist = int(max(random.normalvariate(pair_dist, pair_sd), 5))
        # Chose which haplotype to use as the ref
        ref_hap = random.randrange(len(haplotypes))
        seq_range = len(haplotypes[ref_hap]) - this_pair_dist + 1
        frag_start = random.randrange(seq_range)
        if isinstance(spos, list):
            spos.append(frag_start)
        frag = haplotypes[ref_hap][frag_start:frag_start + this_pair_dist]
        if random.random() < 0.5:
            frag = revcomp(frag)
        frag = mutate_sequence(frag, rate=errrate)
        if len(frag) >= readlen:
            r1 = frag[:readlen]
            r2 = revcomp(frag[-readlen:])
        else:
            r1 = frag + adaptors[0]
            r1.ljust(readlen, 'N')
            r2 = revcomp(frag) + adaptors[1]
            r1.ljust(readlen, 'N')
        for i, r in enumerate([r1, r2]):
            name = '{:d}-{:d}/{:d}'.format(pair_num, frag_start, i)
            yield fq(name, r, qual)

def wgsim(genome, num_pairs, file=sys.stdout, mutrate=0.001, indelfrac=0.1,
          indexext=0.3, haploid=False, readlen=101, errrate=0.01,
          pair_dist=500, pair_sd=30, seed=None):
    haps = mutate_ref(genome, mutrate, indelfrac, indelext, haploid, seed)
    for read in sim_reads(haps, num_pairs, errrate, pair_dist, pair_sd, seed):
        file.write(read.encode('ascii'))
