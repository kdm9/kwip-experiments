# Mapping of population name: divergence
GENOMES = ['A', 'B', 'C', 'D']
SAMPLES = [str(i) for i in range(1, 5)]
SAMPLE_SEEDS = {g: {s: str(i + len(SAMPLES) * j + 1) for j, s in enumerate(SAMPLES)}
                for i, g in enumerate(GENOMES)}
GENOME_SIZE = int(1e7)
COVERAGES = [1, 2, 5, 10, 20, 50]
SCALES = ['0.001', '0.01', '0.1']
READ_NUMS = {cov: int(GENOME_SIZE * cov / 200) for cov in COVERAGES}
HASH_SIZE = "1e9"
METRICS = ['wip', 'ip']


rule all:
    input:
        expand("data/genomes-{s}/{g}.fasta", g=GENOMES, s=SCALES),
        expand("data/samples/{genome}-{scale}-{sample}_{cov}x_il.fastq.gz",
               genome=GENOMES, scale=SCALES, sample=SAMPLES, cov=COVERAGES),
        expand("data/kwip/{cov}x-{scale}.stat", cov=COVERAGES, scale=SCALES),
        expand("data/kwip/{cov}x-{scale}-{metric}.{ext}", cov=COVERAGES,
               scale=SCALES, metric=METRICS, ext=["dist", "kern"]),

rule all_genomes:
    input:
        "tree.nwk"
    output:
        temp("data/merged_genomes-{scale}.fasta")
    params:
        length=str(GENOME_SIZE),
        seed="23",
        scale=lambda w: w.scale,
    log:
        "data/log/seqgen-{scale}.log"
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
        "data/merged_genomes-{scale}.fasta"
    output:
        expand("data/genomes-{{scale}}/{genome}.fasta", genome=GENOMES),
    params:
        dir=lambda w: "data/genomes-" + w.scale
    log:
        "data/log/splitfa-{scale}.log"
    shell:
        "./splitfa.py"
        " {input}"
        " {params.dir}/"
        " >{log} 2>&1"


rule alignment:
    input:
        expand("data/genomes-{{scale}}/{genome}.fasta", genome=GENOMES)
    output:
        "data/aligned_genomes-{scale}.fasta"
    log:
        "data/log/alignment-{scale}.log"
    shell:
        "cat {input}"
        " | mafft  --nuc -"
        " >{output}"
        " 2>{log}"


rule samples:
    input:
        "data/genomes-{scale}/{genome}.fasta",
    output:
        r1=temp("/dev/shm/{genome}-{scale}-{sample}_{cov}x_R1.fastq"),
        r2=temp("/dev/shm/{genome}-{scale}-{sample}_{cov}x_R2.fastq")
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
        r1="/dev/shm/{genome}-{scale}-{sample}_{cov}x_R1.fastq",
        r2="/dev/shm/{genome}-{scale}-{sample}_{cov}x_R2.fastq"
    priority:
        10
    output:
        "data/samples/{genome}-{scale}-{sample}_{cov}x_il.fastq.gz"
    log:
        "data/log/join/{genome}-{scale}-{sample}_{cov}x.log"
    shell:
        "pairs join"
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


rule kwip:
    input:
        sorted(expand("data/hashes/{genome}-{{scale}}-{sample}_{{cov}}x.ct.gz",
               genome=GENOMES, sample=SAMPLES))
    output:
        d="data/kwip/{cov}x-{scale}-{metric}.dist",
        k="data/kwip/{cov}x-{scale}-{metric}.kern"
    params:
        metric= lambda w: '-U' if w.metric == 'ip' else ''
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
        expand("data/hashes/{genome}-{{scale}}-{sample}_{{cov}}x.ct.gz",
               genome=GENOMES, sample=SAMPLES)
    output:
        "data/kwip/{cov}x-{scale}.stat"
    log:
        "data/log/kwip-entvec/{cov}x-{scale}.log"
    threads:
        24
    shell:
        "kwip-entvec"
        " -o {output}"
        " -t {threads}"
        " {input}"
        " >{log} 2>&1"
