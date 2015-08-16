from __future__ import print_function
import random
import sys
import textwrap


NTs = ['A', 'C', 'G', 'T']


def make_rand_genome(mbp=10, n_freq=0.00001):
    gen = []
    bp = int(mbp * 1000000)
    while len(gen) < bp:
        if random.random() < n_freq:
            gen.append('N')
        else:
            gen.append(random.choice(NTs))
    return ''.join(gen)


def seq_to_fasta_rec(seq, name, wrap=False):
    lines = []
    lines.append('>{}'.format(name))
    if wrap:
        lines.extend(textwrap.wrap(seq, width=80))
    else:
        lines.append(seq)
    return '\n'.join(lines)


def mutate_sequence(seq, rate=0.001):
    nts = list(seq)
    for idx, nt in enumerate(nts):
        if random.random() < rate:
            new_nt = nt
            while new_nt == nt:
                new_nt = random.choice(NTs)
            nts[idx] = new_nt
    return ''.join(nts)


def biforcating_sequences(seq, levels=4, av_rate=0.001,
                          sd_rate=0.00001):
    a_rate = max(random.normalvariate(av_rate, sd_rate), 0.0)
    b_rate = max(random.normalvariate(av_rate, sd_rate), 0.0)
    a = mutate_sequence(seq, a_rate)
    b = mutate_sequence(seq, b_rate)
    print("at level", levels, "(a rate:", a_rate, "b rate:", b_rate, ")",
          file=sys.stderr)
    if levels > 0:
        return (biforcating_sequences(a, levels - 1, av_rate, sd_rate),
                biforcating_sequences(b, levels - 1, av_rate, sd_rate))
    else:
        return (a, b)


def flatten(lst):
    for item in lst:
        if isinstance(item, list) or isinstance(item, tuple):
            for x in flatten(item):
                yield x
        else:
            yield item


def print_multifasta(seqlist, file=sys.stdout):
    for i, seq in enumerate(seqlist):
        print(seq_to_fasta_rec(seq, str(i + 1)), file=file)
