#!/bin/bash
mkdir -p clogs
mkdir -p jobs

set -xe

HERE=$PWD
#for cfg in sets/rice/*.json
for cfg in sets/rice/0.json
do
	job=$(basename $cfg .json)
	jobscript=jobs/${job}.job
	cat >${jobscript} <<EOF
#PBS -l ncpus=16
#PBS -l walltime=24:00:00
#PBS -P xe2
#PBS -q normal
#PBS -l other=gdata1
#PBS -l mem=126GB
#PBS -l wd
#PBS -o ./clogs
#PBS -e ./clogs

source ~/.profile
source /g/data1/xe2/.profile
set -xe

workdir=/short/xe2/kwip-experiments/realdata/runs/$job
mkdir -p $workdir
cd $workdir
mkdir -p data
mkdir -p logs
ln -s /g/data1/xe2/datasets/3000-rice/sra/ data/sra

snakemake -j 16 --configfile $HERE/$cfg --snakefile $HERE/Snakefile >logs/${job}.log 2>&1
EOF
	qsub $jobscript
done
