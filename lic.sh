set -xe
ls data/*.fq.gz |parallel --verbose \
	load-into-counting.py -x 5e8 -N 1 -k 20 -b -s tsv {.}.ct.gz {}
