# Mapping of population name: divergence
import string

Ne = 2048   # Effective population size
Npop = 8    # populations
Nloci = 5   # Chromosomes
N = 48      # number of samples from population
REPS = 2    # Replicate runs per sample


# Sample names
SAMPLES = [string.ascii_uppercase[i] for i in range(1, 5)]
COVERAGES = [1, 5, 10, 20,]

GENOME_SIZE = int(1e6)
HASH_SIZE = "5e7"
METRICS = ['wip', 'ip']


rule all:
    input:
        expand("data/pops/P{psz}-p{np}.nwk", psz=POPSIZES, np=NPOP),
        expand("data/genomes/N{N}-P{nhap}-p{npop}/", N=Ns, nhap=POPSIZES, npop=NPOP),

rule clean:
    shell:
        "rm -rf data"


def scrm_args(wc):
    migration = 0.3
    nhap = int(wc.nhap)
    nrep = 10 #int(wc.nrep)
    npop = int(wc.npop)
    args = "{:d} {:d} ".format(nhap, nrep)
    args += "-T "

    perpop = int(nhap / npop)
    pops = " ".join([str(perpop) for _ in range(npop)])
    args += " -I {:d} {:s} {} ".format(npop, pops, migration)

    return args

rule population:
    output:
        "data/pops/P{nhap}-p{npop}.nwk"
    params:
        args=scrm_args
    log:
        "data/log/scrm-P{nhap}-p{npop}.log"
    shell:
        "scrm"
        " {params.args}"
        " 2>{log}"
        "| grep '(' "  # Grep for (, which are only in a newick tree output line
        " >{output}"


rule sampled_tree:
    input:
        "data/pops/P{nhap}-p{npop}.nwk"
    output:
        "data/samples/N{N}-P{nhap}-p{npop}.nwk"
    params:
        N=lambda w: w.N
    log:
        "data/log/scrm-P{nhap}-p{npop}.log"
    shell:
        "python ./random_subtree.py"
        " {input}"
        " {params.N}"
        " 2>{log}"
        " >{output}"


rule all_genomes:
    input:
        "data/samples/N{N}-P{nhap}-p{npop}.nwk"
    output:
        temp("data/merged_genomes-N{N}-P{nhap}-p{npop}.fasta")
    params:
        length=str(GENOME_SIZE),
        seed="23",
        scale=SCALE,
    log:
        "data/log/seqgen.log"
    shell:
        "seq-gen"
        " -mF84"
        " -s{params.scale}"
        " -l{params.length}"
        " -z{params.seed}"
        " -oF"
        " {input}"
        " >{output} 2>{log}"


rule genomes:
    input:
        "data/merged_genomes-N{N}-P{nhap}-p{npop}.fasta"
    output:
        "data/genomes/N{{N}}-P{{nhap}}-p{{npop}}/*.fasta",
        #expand("data/genomes/N{{N}}-P{{nhap}}-p{{npop}}/{genome}.fasta", genome=GENOMES),
    params:
        dir=lambda w: "data/genomes/N{}-P{}-p{}".format(w.N, w.nhap, w.npop)
    log:
        "data/log/splitfa-N{N}-P{nhap}-p{npop}.log"
    shell:
        "./splitfa.py"
        " {input}"
        " {params.dir}/"
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
