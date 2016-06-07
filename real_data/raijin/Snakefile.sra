SETS = {}
ALLRUNS = set()
RUNPROJECTS = {}
for project, setruns in config.items():
    for setname, runs in setruns.items():
        for r in runs:
            ALLRUNS.add(r)
ALLRUNS = list(ALLRUNS)
PROJECT = []
for run in ALLRUNS:
    PROJECT.append(config['runprojects']['run'])


## BEGIN RULES
rule all:
    input:
        expand("data/sra/{project}/{run}.sra", zip, PROJECT, ALLRUNS)

rule sra:
    output:
        "data/sra/{project}/{run}.sra",
    log:
        "data/log/getrun/{run}.log"
    params:
        srr=lambda w: w.run
    shell:
        "get-run.py"
        "   -d data/sra"
        "   -s "
        "   -i {wildcards.run}"
        " >{log} 2>&1"
