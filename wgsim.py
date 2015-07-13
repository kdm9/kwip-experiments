import random

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


def sim_reads(haplotypes, num_pairs, readlen=101, errrate=0.01, pair_dist=500,
              pair_sd=30, seed=None):
    random.seed(seed)

    qual_char = 'I'  # Maximum score in PHRED 33, for errrate == 0
    if errrate > 0:
        qual_char = chr(int(-10.0 * log10(errrate)) + 33)

    for pair_num in range(num_pairs):
        pass


def wgsim(genome, num_reads, readlen=101, errrate=0.01, mutrate=0.001,
          indelfrac=0.1, indexext=0.3, pair_dist=500, pair_sd=30, num_hap=2,
          seed=None):
    pass
