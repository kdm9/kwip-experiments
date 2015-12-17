# Mapping of population name: divergence
import string
import itertools
from math import log, ceil
from snakemake.utils import format as snakefmt
import random

#Ne = 2048   # Effective population size
Ne = 512   # Effective population size
Npop = 8    # populations
Nloci = 1   # Chromosomes
#N = 48      # number of samples from population
N = 12      # number of samples from population
REPS = 3    # Replicate runs per sample
M = 0.3     # Migration Rate


# Sample names
labels = [
    "".join(x)
    for x in itertools.product(string.ascii_uppercase, repeat=int(ceil(log(N, 26))))
]
GENOMES = [labels[i] for i in range(N)]
SAMPLES = list(range(REPS))
COVERAGES = [1, 5, 10, 20,]
COVERAGES = [5, 20,]
READ_NUMS = {}
GENOME_SIZE = int(1e6)
HASH_SIZE = "5e7"
SCALE = 0.01  # per base subs/indel rate


random.seed(234)
for g in GENOMES:
    READ_NUMS[g] = {}
    for s in SAMPLES:
        READ_NUMS[g][str(s)] = {}
        for c in COVERAGES:
            cov = random.gauss(c, c* 0.3)
            nread = (cov * GENOME_SIZE) / 200
            READ_NUMS[g][str(s)][str(c)] = int(nread)

METRICS = ['wip', 'ip']

KWIPS = {
    '0.1.8': './bin/kwip-0.1.8',
    '0.1.10': './bin/kwip-0.1.10',
    'l1norm': './bin/kwip-l1',
}


rule all:
    input:
        expand("data/kwip/{kwip}/{kwip}-{cov}x-{metric}.{ext}", cov=COVERAGES,
                metric=METRICS, ext=['dist', 'kern'], kwip=KWIPS),
        expand("data/kwip/{cov}x.stat", cov=COVERAGES),
        #expand("data/kwip/{kwip}/{cov}x.stat", cov=COVERAGES, kwip=KWIPS),

rule clean:
    shell:
        "rm -rf data .snakemake"


def scrm_args(wc):
    args = "{:d} {:d} ".format(Ne, Nloci)
    args += "-T "
    perpop = int(Nloci / Npop)
    pops = " ".join([str(perpop) for _ in range(Npop)])
    args += " -I {:d} {:s} {} ".format(Npop, pops, M)
    return args

rule population:
    output:
        "data/population.nwk"
    params:
        args=scrm_args
    log:
        "data/log/scrm.log"
    shell:
        "scrm"
        " {params.args}"
        " 2>{log}"
        "| grep '(' "  # Grep for (, which are only in a newick tree output line
        " >{output}"


rule dawgtree:
    input:
        "data/population.nwk"
    output:
        "data/sample.nwk"
    params:
        N=str(N)
    shell:
        "scripts/random_subtree.py"
        " -n {params.N}"
        " -o {output}"
        " {input}"


rule paramfile:
    input:
        "data/sample.nwk"
    output:
        "data/dawg.params"
    params:
        length=str(GENOME_SIZE),
    run:
        with open(input[0])  as fh:
            tree = fh.read().strip()
        dawg = """
        [Tree]
        Tree = "{tree}"
        Scale = {scale:.20f}

        [Indel]
        Model.Ins = Geo
        Model.Del = Geo
        Rate.Ins = 0.005
        Rate.Del = 0.005

        [Subst]
        Model = F84
        Freqs  = 0.3, 0.2, 0.2, 0.3
        Params = 2.5

        [Root]
        Length = {length}
        """.format(length=str(GENOME_SIZE), tree=tree, scale=SCALE * 1e-6)
        with open(output[0], 'w') as fh:
            print(dawg, file=fh)

rule all_genomes:
    input:
        "data/dawg.params",
    output:
        "data/all_genomes.fasta",
    log:
        "data/log/dawg.log"
    shell:
        "dawg2"
        " {input}"
        " -o {output}"
        " >{log} 2>&1"


rule genomes:
    input:
        "data/all_genomes.fasta"
    output:
        expand("data/genomes/{g}.fasta", g=GENOMES)
    params:
        dir="data/genomes/"
    log:
        "data/log/splitfa.log"
    shell:
        "scripts/splitfa.py"
        " {input}"
        " {params.dir}"
        " >{log} 2>&1"


rule samples:
    input:
        "data/genomes/{genome}.fasta",
    output:
        r1=temp("data/tmp/{genome}-{sample}_{cov}x_R1.fastq"),
        r2=temp("data/tmp/{genome}-{sample}_{cov}x_R2.fastq")
    params:
        seed=lambda w: str(int(hash(w.genome + w.sample) % 999983)),
        rn=lambda w: str(READ_NUMS[w.genome][w.sample][w.cov])
    log:
        "data/log/samples/{genome}-{sample}_{cov}x.log"
    shell:
        "mason_simulator"
        " -ir {input}"
        #" --illumina-prob-mismatch-scale 2.0"
        " --illumina-read-length 101"
        " -o {output.r1}"
        " -or {output.r2}"
        " --seed {params.seed}"
        " -n {params.rn}"
        " >{log} 2>&1"


rule ilfq:
    input:
        r1="data/tmp/{genome}-{sample}_{cov}x_R1.fastq",
        r2="data/tmp/{genome}-{sample}_{cov}x_R2.fastq"
    priority:
        10
    output:
        "data/samples/{genome}-{sample}_{cov}x_il.fastq.gz"
    log:
        "data/log/join/{genome}-{sample}_{cov}x.log"
    shell:
        "interleave-reads.py"
        " {input.r1}"
        " {input.r2}"
        " 2>>{log}"
        " | ./bin/trimit"
        " -q 28"
        " 2>/dev/null"
        " | gzip > {output}"
        " 2>>{log}"


rule hash:
    input:
        "data/samples/{genome}-{sample}_{cov}x_il.fastq.gz"
    output:
        "data/hashes/{genome}-{sample}_{cov}x.ct.gz"
    params:
        x=HASH_SIZE,
        N='1',
        k='20',
    log:
        "data/log/hashes/{genome}-{sample}_{cov}x.log"
    shell:
        "load-into-counting.py"
        " -N {params.N}"
        " -x {params.x}"
        " -k {params.k}"
        " -b"
        " -s tsv"
        " {output}"
        " {input}"
        " >{log} 2>&1"


rule kwip:
    input:
        sorted(expand("data/hashes/{genome}-{sample}_{{cov}}x.ct.gz",
               genome=GENOMES, sample=SAMPLES))
    output:
        d="data/kwip/{kwip}/{kwip}-{cov}x-{metric}.dist",
        k="data/kwip/{kwip}/{kwip}-{cov}x-{metric}.kern"
    params:
        metric= lambda w: '-U' if w.metric == 'ip' else '',
        kwip=lambda w: KWIPS[w.kwip],
    log:
        "data/log/kwip/{kwip}/{cov}x-{metric}.log"
    threads:
        24
    shell:
        "{params.kwip}"
        " {params.metric}"
        " -d {output.d}"
        " -k {output.k}"
        " -t {threads}"
        " {input}"
        " >{log} 2>&1"


rule kwip_stats:
    input:
        expand("data/hashes/{genome}-{sample}_{{cov}}x.ct.gz",
               genome=GENOMES, sample=SAMPLES)
    output:
        "data/kwip/{cov}x.stat"
    log:
        "data/log/kwip-entvec/{cov}x.log"
    threads:
        24
    shell:
        "kwip-entvec"
        " -o {output}"
        " -t {threads}"
        " {input}"
        " >{log} 2>&1"


# rule kwip_stats:
#     input:
#         expand("data/hashes/{genome}-{sample}_{{cov}}x.ct.gz",
#                genome=GENOMES, sample=SAMPLES)
#     output:
#         "data/kwip/{kwip}/{cov}x.stat"
#     params:
#         kwip=lambda w: KWIPS[w.kwip],
#     log:
#         "data/log/kwip-entvec/{kwip}/{cov}x.log"
#     threads:
#         24
#     shell:
#         "{params.kwip}"
#         " -o {output}"
#         " -t {threads}"
#         " {input}"
#         " >{log} 2>&1"
