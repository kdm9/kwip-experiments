NP=8
NH=4096
NR=10
PP=$(($NH / $NP))
POPS=""
N=48

for x in `seq 1 $NP`
do
    POPS="$POPS $PP"
done

scrm $NH $NR -T -I $NP $POPS 0.3 | grep '(' > data/popn_tree.nwk

python ./random_subtree_cog.py -n $N data/popn_tree.nwk >data/sample_tree.nwk

seq-gen -t3.5 -mHKY -f 0.2 0.3 0.3 0.2 -s0.000001 -l1000000 data/sample_tree.nwk >data/genomes.fa
