# Mapping of population name: divergence
import string
import itertools
from math import log, ceil
from snakemake.utils import format as snakefmt

Ne = 2048   # Effective population size
Npop = 8    # populations
Nloci = 1   # Chromosomes
N = 48      # number of samples from population
REPS = 2    # Replicate runs per sample
M = 0.3     # Migration Rate


# Sample names
labels = [
    "".join(x)
    for x in itertools.product(string.ascii_uppercase, repeat=int(ceil(log(N, 26))))
]
SAMPLES = [labels[i] for i in range(N)]
COVERAGES = [1, 5, 10, 20,]

GENOME_SIZE = int(5e5)
HASH_SIZE = "5e7"
METRICS = ['wip', 'ip']


rule all:
    input:
        expand("data/genomes/{g}.fasta", g=SAMPLES),

rule clean:
    shell:
        "rm -rf data"


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
        "data/dawg.tree"
    params:
        N=str(N)
    shell:
        "scripts/random_subtree_cog.py"
        " -n {params.N}"
        " {input}"

rule paramfile:
    output:
        "data/dawg.params"
    params:
        length=str(GENOME_SIZE),
    run:
        dawg = """
        Length = {{{length}, }}
        Seed = 23
        TreeScale = 2e-5
        Model = 'F84'
        Params = {{4,}}
        Lambda = 0.1
        File = 'data/all_genomes.fasta'
        Format = 'Fasta'
        GapSingleChar = true
        GapModel = 'NB'
        GapParams = {{1, 0.5}}
        """.format(length=", ".join(map(str, [GENOME_SIZE] * Nloci)))
        with open(output[0], 'w') as fh:
            print(dawg, file=fh)

rule all_genomes:
    input:
        "data/dawg.params",
        "data/dawg.tree",
    output:
        temp("data/all_genomes.fasta")
    log:
        "data/log/dawg.log"
    shell:
        "dawg"
        " -c"
        " {input}"
        " >{log} 2>&1"


rule genomes:
    input:
        "data/all_genomes.fasta"
    output:
        expand("data/genomes/{g}.fasta", g=SAMPLES)
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
        "data/genomes-{scale}/{genome}.fasta",
    output:
        r1=temp("data/tmp/{genome}-{scale}-{sample}_{cov}x_R1.fastq"),
        r2=temp("data/tmp/{genome}-{scale}-{sample}_{cov}x_R2.fastq")
    params:
        seed=lambda w: SAMPLE_SEEDS[w.genome][w.sample],
        rn=lambda w: str(READ_NUMS[int(w.cov)])
    log:
        "data/log/samples/{genome}-{scale}-{sample}_{cov}x.log"
    shell:
        "mason_simulator"
        " -ir {input}"
        " --illumina-prob-mismatch-scale 2.0"
        " --illumina-read-length 101"
        " -o {output.r1}"
        " -or {output.r2}"
        " --seed {params.seed}"
        " -n {params.rn}"
        " >{log} 2>&1"


rule ilfq:
    input:
        r1="data/tmp/{genome}-{scale}-{sample}_{cov}x_R1.fastq",
        r2="data/tmp/{genome}-{scale}-{sample}_{cov}x_R2.fastq"
    priority:
        10
    output:
        "data/samples/{genome}-{scale}-{sample}_{cov}x_il.fastq.gz"
    log:
        "data/log/join/{genome}-{scale}-{sample}_{cov}x.log"
    shell:
        "interleave-reads.py"
        " {input.r1}"
        " {input.r2}"
        " 2>>{log}"
        " | sickle se "
        " -f /dev/stdin"
        " -t sanger"
        " -o /dev/stdout"
        " -n"
        " 2>>{log}"
        " | gzip > {output}"
        " 2>>{log}"


rule hash:
    input:
        "data/samples/{genome}-{scale}-{sample}_{cov}x_il.fastq.gz"
    output:
        "data/hashes/{genome}-{scale}-{sample}_{cov}x.ct.gz"
    params:
        x=HASH_SIZE,
        N='1',
        k='20',
    log:
        "data/log/hashes/{genome}-{scale}-{sample}_{cov}x.log"
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


# rule kwip:
#     input:
#         sorted(expand("data/hashes/{genome}-{{scale}}-{sample}_{{cov}}x.ct.gz",
#                genome=GENOMES, sample=SAMPLES))
#     output:
#         d="data/kwip/{cov}x-{scale}-{metric}.dist",
#         k="data/kwip/{cov}x-{scale}-{metric}.kern"
#     params:
#         metric= lambda w: '-U' if w.metric == 'ip' else ''
#     log:
#         "data/log/kwip/{cov}x-{scale}-{metric}.log"
#     threads:
#         24
#     shell:
#         "kwip"
#         " {params.metric}"
#         " -d {output.d}"
#         " -k {output.k}"
#         " -t {threads}"
#         " {input}"
#         " >{log} 2>&1"
# 
# 
# rule kwip_stats:
#     input:
#         expand("data/hashes/{genome}-{{scale}}-{sample}_{{cov}}x.ct.gz",
#                genome=GENOMES, sample=SAMPLES)
#     output:
#         "data/kwip/{cov}x-{scale}.stat"
#     log:
#         "data/log/kwip-entvec/{cov}x-{scale}.log"
#     threads:
#         24
#     shell:
#         "kwip-entvec"
#         " -o {output}"
#         " -t {threads}"
#         " {input}"
#         " >{log} 2>&1"
