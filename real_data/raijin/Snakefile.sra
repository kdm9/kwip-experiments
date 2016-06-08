SETS = {}
ALLRUNS = set()
RUNPROJECTS = {}
for project, setruns in config.items():
    for setname, runs in setruns.items():
        for r in runs:
            ALLRUNS.add(r)
            RUNPROJECTS[r] = project

## BEGIN RULES

rule all:
    input:
        ["data/sra/{project}/{run}.sra".format(run=run, project=project) for run, project in RUNPROJECTS.items()]

rule sra:
    output:
        "data/sra/{project}/{run}.sra",
    log:
        "data/log/getrun/{project}-{run}.log"
    params:
    shell:
        "get-run.py"
        "   -d data/sra/{wildcards.project}"
        "   -s "
        "   -i {wildcards.run}"
        " >{log} 2>&1"
