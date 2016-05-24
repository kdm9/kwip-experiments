#PBS -l ncpus=1
#PBS -l walltime=24:00:00
#PBS -P xe2
#PBS -q express
#PBS -l other=gdata1
#PBS -l mem=10GB
#PBS -l wd
#PBS -o ./clogs
#PBS -e ./clogs

mkdir -p clogs/

( \
snakemake 					\
	-j 200					\
	--configfile sets/rice/all.json		\
       	--cluster-config cluster.json		\
       	--js jobscript.sh 			\
	--cluster 'qsub -q {cluster.queue} -l ncpus={threads} -l walltime={cluster.time} -l mem={cluster.mem} {cluster.miscargs}' \
) > snakemake.log 2>&1

