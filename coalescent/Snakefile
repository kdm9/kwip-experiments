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


SEED = 435
random.seed(SEED)

# Sample names
labels = [
    "".join(x)
    for x in itertools.product(string.ascii_uppercase, repeat=int(ceil(log(N, 26))))
]

GENOMES = [labels[i] for i in range(N)]
SAMPLES = list(range(REPS))
COV_VAR = [
    (0.1, 0.01),
    (0.5, 0.01),
    (1, 0.01),
    (2, 0.01),
    (5, 0.01),
    (10, 0.01),
    (15, 0.01),
    (30, 0.01),
    (50, 0.01),
    (100, 0.01),
    (30, 0.05),
    (30, 0.001),
]

COVERAGE_CV = 0.3
GENOME_SIZE = int(1e7)
HASH_SIZE = "5e8"
SCALE = 0.01  # per base subs/indel rate


# All sim params
READ_NUMS = {}

for g in GENOMES:
    READ_NUMS[g] = {}
    for s in SAMPLES:
        READ_NUMS[g][str(s)] = {}
        for c, v in COV_VAR:
            cov = random.gauss(c, c * COVERAGE_CV)
            nread = (cov * GENOME_SIZE) / 200
            try:
                READ_NUMS[g][str(s)][str(c)][str(v)] = int(nread)
            except KeyError:
                READ_NUMS[g][str(s)][str(c)] = {str(v): int(nread)}


METRICS = ['wip', 'ip']


rule all:
    input:
        expand(expand("data/hashes/{cov}x-{scale}/{{genome}}-{{sample}}.ct.gz", zip,
                      cov=map(lambda x: x[0], COV_VAR),
                      scale=map(lambda x: x[1], COV_VAR)),
               genome=GENOMES, sample=SAMPLES),
        expand(expand("data/kwip/{cov}x-{scale}-{{metric}}.{{ext}}", zip,
                      cov=map(lambda x: x[0], COV_VAR),
                      scale=map(lambda x: x[1], COV_VAR)),
               metric=METRICS, ext=['dist', 'kern']),
        expand("data/kwip/{cov}x-{scale}.stat", zip,
               cov=map(lambda x: x[0], COV_VAR),
               scale=map(lambda x: x[1], COV_VAR)),


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
        "data/dawg-{scale}.params",
    params:
        length=str(GENOME_SIZE),
        scale=lambda wc: "{:.20f}".format(float(wc.scale) * 1e-6)
    run:
        with open(input[0])  as fh:
            tree = fh.read().strip()
        dawg = """
        [Tree]
        Tree = "{tree}"
        Scale = {scale}

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
        """.format(length=str(GENOME_SIZE), tree=tree, scale=params.scale)
        with open(output[0], 'w') as fh:
            print(dawg, file=fh)


rule all_genomes:
    input:
        "data/dawg-{scale}.params",
    output:
        "data/all_genomes-{scale}.fasta"
    log:
        "data/log/dawg-{scale}.log"
    shell:
        "dawg2"
        " {input}"
        " -o {output}"
        " >{log} 2>&1"


rule genomes:
    input:
        "data/all_genomes-{scale}.fasta"
    output:
        expand("data/genomes/{{scale}}/{g}.fasta", g=GENOMES)
    params:
        dir="data/genomes/{scale}/"
    log:
        "data/log/splitfa-{scale}.log"
    shell:
        "scripts/splitfa.py"
        " {input}"
        " {params.dir}"
        " >{log} 2>&1"


rule samples:
    input:
        "data/genomes/{scale}/{genome}.fasta",
    output:
        r1=temp("data/temp/{cov}x-{scale}/{genome}-{sample}_R1.fastq"),
        r2=temp("data/temp/{cov}x-{scale}/{genome}-{sample}_R2.fastq"),
    params:
        seed=lambda w: str(int(hash(w.genome + w.sample) % 999983)),
        rn=lambda w: str(READ_NUMS[w.genome][w.sample][w.cov][w.scale])
    log:
        "data/log/samples/{cov}x-{scale}/{genome}-{sample}.log"
    shell:
        "mason_simulator"
        " -ir {input}"
        " --illumina-read-length 101"
        " -o {output.r1}"
        " -or {output.r2}"
        " --seed {params.seed}"
        " -n {params.rn}"
        " >{log} 2>&1"


rule ilfq:
    input:
        r1="data/temp/{cov}x-{scale}/{genome}-{sample}_R1.fastq",
        r2="data/temp/{cov}x-{scale}/{genome}-{sample}_R2.fastq"
    priority:
        10
    output:
        temp("data/temp/{cov}x-{scale}/{genome}-{sample}_il.fastq")
    log:
        "data/log/join/{cov}x-{scale}/{genome}-{sample}.log"
    shell:
        "interleave-reads.py"
        " {input.r1}"
        " {input.r2}"
        " 2>>{log}"
        " | trimit"
        " -q 28"
        " 2>/dev/null"
        " > {output}"


rule hash:
    input:
        "data/temp/{cov}x-{scale}/{genome}-{sample}_il.fastq"
    output:
        "data/hashes/{cov}x-{scale}/{genome}-{sample}.ct.gz"
    priority:
        11
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
        expand("data/hashes/{{cov}}x-{{scale}}/{genome}-{sample}.ct.gz",
               genome=GENOMES, sample=SAMPLES),
    output:
        d="data/kwip/{cov}x-{scale}-{metric}.dist",
        k="data/kwip/{cov}x-{scale}-{metric}.kern"
    params:
        metric= lambda w: '-U' if w.metric == 'ip' else '',
    log:
        "data/log/kwip/{cov}x-{scale}-{metric}.log"
    threads:
        24
    shell:
        "kwip"
        " {params.metric}"
        " -d {output.d}"
        " -k {output.k}"
        " -t {threads}"
        " {input}"
        " >{log} 2>&1"


rule kwip_stats:
    input:
        expand("data/hashes/{{cov}}x-{{scale}}/{genome}-{sample}.ct.gz",
               genome=GENOMES, sample=SAMPLES),
    output:
        "data/kwip/{cov}x-{scale}.stat"
    log:
        "data/log/kwip-stats/{cov}x-{scale}.log"
    threads:
        24
    shell:
        "kwip-stats"
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
