set -xe
bwa index data/genome.fa
bwa mem -p data/genome.fa data/reads.fq.gz | \
	samtools view -Su - | \
	samtools sort -o data/reads_s.bam -T ~/tmp

samtools index data/reads_s.bam
