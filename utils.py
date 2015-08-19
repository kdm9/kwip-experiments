from __future__ import print_function
import random
import sys
import textwrap
import os
import subprocess
import skbio
import gzip

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
    if levels > 1:
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


def run_cmd(cmd, quiet=True):
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in p.stdout:
        if quiet:
            continue
        print(line.decode('utf-8'), end='')
    return p.wait()

def wgsim(N, ref, ilfq, rate=0.0000001, err_rate=0.001, quiet=True):
    if N <= 10000000:
        tmp='/dev/shm'
    else:
        tmp='/tmp'
    rand = random.randint(0, 100000000000)
    r1_file = '{tmp}/{r}_r1.fq'.format(tmp=tmp, r=rand)
    r2_file = '{tmp}/{r}_r2.fq'.format(tmp=tmp, r=rand)
    try:
        wgs = "wgsim -S {rand} -e {err} -1 101 -2 101 -r {rate:.12f} -N {N} {ref} {r1} {r2}".format(
                rand=rand, err=err_rate, N=N, ref=ref, rate=rate, r1=r1_file, r2=r2_file)
        if not quiet:
            print(wgs)
        run_cmd(wgs)
        pairs = 'pairs join {r1} {r2} | pigz > {ilfq}'.format(r1=r1_file, r2=r2_file, ilfq=ilfq)
        run_cmd(pairs, False)
    finally:
        # nuke the temp files
        try:
            os.remove(r1_file)
            os.remove(r2_file)
        except:
            pass

