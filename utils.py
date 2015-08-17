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



def wgsim(N, ref, ilfq, rate=0.0000001, quiet=True):
    if N <= 10000000:
        tmp='/dev/shm'
    else:
        tmp='/tmp'
    rand = random.randint(0, 100000000000)
    r1_file = '{tmp}/{r}_r1.fq'.format(tmp=tmp, r=rand)
    r2_file = '{tmp}/{r}_r2.fq'.format(tmp=tmp, r=rand)
    try:
        wgs = "wgsim -S {rand} -1 101 -2 101 -r {rate:.12f} -N {N} {ref} {r1} {r2}".format(
                rand=rand, N=N, ref=ref, rate=rate, r1=r1_file, r2=r2_file)
        if not quiet:
            print(wgs)
        p = subprocess.Popen(wgs, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout:
            if quiet:
                continue
            print(line.decode('ascii'), end='')
        p.wait()
        with open(ilfq, 'w') as il_fh:
            r1 = []
            r2 = []
            with open(r1_file) as r1_fh, open(r2_file) as r2_fh:
                for i, (l1, l2) in enumerate(zip(r1_fh, r2_fh)):
                    if i % 4 == 0:
                        for line in r1 + r2:
                            il_fh.write(line)
                        r1 = []
                        r2 = []
                    r1.append(l1)
                    r2.append(l2)
            #r1_fh = skbio.read(r1_file, format='fastq', variant='sanger')
            #r2_fh = skbio.read(r2_file, format='fastq', variant='sanger')
            #for r1, r2 in zip(r1_fh, r2_fh):
            #    skbio.write(r1, format='fastq', variant='sanger', into=il_fh)
            #    skbio.write(r2, format='fastq', variant='sanger', into=il_fh)
    finally:
        # nuke the temp files
        try:
            os.remove(r1_file)
            os.remove(r2_file)
        except:
            pass