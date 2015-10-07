# Mapping of population name: divergence
GENOMES = ['A', 'B', 'C', 'D']
SAMPLES = [str(i) for i in range(1, 5)]
SAMPLE_SEEDS = {g: {s: str(i + len(SAMPLES) * j + 1) for j, s in enumerate(SAMPLES)}
                for i, g in enumerate(GENOMES)}
GENOME_SIZE = int(1e6)
COVERAGES = [1, 2, 5, 10, 50, 100]
READ_NUMS = {cov: int(GENOME_SIZE * cov / 200) for cov in COVERAGES}
HASH_SIZES = ["1e9",]
METRICS = ['wip', 'ip']


rule all:
    input:
        expand("data/genomes/{g}.fasta", g=GENOMES),
        expand("data/samples/{genome}-{sample}_{cov}x_il.fastq.gz",
               genome=GENOMES, sample=SAMPLES, cov=COVERAGES),
        expand("data/kwip/{cov}x.stat", cov=COVERAGES),
        expand("data/kwip/{cov}x-{metric}.{ext}", cov=COVERAGES,
               metric=METRICS, ext=["dist", "kern"]),
        #"data/aligned_genomes.fasta",


rule all_genomes:
    input:
        "tree.nwk"
    output:
        "data/genomes.fasta"
    params:
        length=str(GENOME_SIZE),
        seed="23"
    log:
        "data/log/seqgen.log"
    shell:
        "seq-gen"
        " -mF84"
        " -s0.001"
        " -l{params.length}"
        " -z{params.seed}"
        " -oF"
        " {input}"
        " >{output} 2>{log}"


rule genomes:
    input:
        "data/genomes.fasta"
    output:
        expand("data/genomes/{genome}.fasta", genome=GENOMES)
    log:
        "data/log/splitgen.log"
    shell:
        "./splitfa.py"
        " {input}"
        " data/genomes/"
        " >{log} 2>&1"


rule alignment:
    input:
        expand("data/genomes/{genome}.fasta", genome=GENOMES)
    output:
        "data/aligned_genomes.fasta"
    log:
        "data/log/alignment.log"
    shell:
        "cat {input}"
        " | mafft  --nuc -"
        " >{output}"
        " 2>{log}"


rule samples:
    input:
        "data/genomes/{genome}.fasta",
    output:
        r1=temp("/dev/shm/{genome}-{sample}_{cov}x_R1.fastq"),
        r2=temp("/dev/shm/{genome}-{sample}_{cov}x_R2.fastq")
    params:
        seed=lambda w: SAMPLE_SEEDS[w.genome][w.sample],
        rn=lambda w: str(READ_NUMS[int(w.cov)])
    log:
        "data/log/samples/{genome}-{sample}_{cov}x.log"
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
        r1="/dev/shm/{genome}-{sample}_{cov}x_R1.fastq",
        r2="/dev/shm/{genome}-{sample}_{cov}x_R2.fastq"
    priority:
        10
    output:
        "data/samples/{genome}-{sample}_{cov}x_il.fastq.gz"
    log:
        "data/log/join/{genome}-{sample}_{cov}x.log"
    shell:
        "pairs join"
        " {input.r1}"
        " {input.r2}"
        " 2>>{log}"
        " | gzip > {output}"
        " 2>>{log}"


rule hash:
    input:
        "data/samples/{genome}-{sample}_{cov}x_il.fastq.gz"
    output:
        "data/hashes/{genome}-{sample}_{cov}x.ct.gz"
    params:
        x='1e8',
        N='1',
        k='20'
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
        d="data/kwip/{cov}x-{metric}.dist",
        k="data/kwip/{cov}x-{metric}.kern"
    params:
        metric= lambda w: '-U' if w.metric == 'ip' else ''
    log:
        "data/log/kwip/{cov}x-{metric}.log"
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
