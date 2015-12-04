================================
Coalescent population simulation
================================

Aim
^^^

This simulation aims to create sample genomes that come from a realistic
simulation of a population. We then simulate whole-genome shotgun sequencing
reads from these genomes at various average coverages and estimate similarity
from the resulting reads with `kWIP`.


Methods
^^^^^^^

We simulate a population using [`scrm`](https://github.com/scrm/scrm). The
simulation we do is of the overall population's history and dynamics, not
simply the sample. The Ne (effective population size) is in the order of a few
thousand (check current Snakefile for detail).

Once we have simulated a population, we randomly sub-sample this tree to create
a population sample. This is designed to simulate scattered geographic sampling
of some natural population. The sample genomes are then generated using Dawg,
using the F84 model of sequence evolution, allowing for an unequal number of
transitions and transversions. Dawg also simulates the occurrence of short
indel variants. We currently simulate a genome size of 1Mbp, to reduce
runtime of the simulation. For the paper we'd probably want to increase that.

The resulting alignment is split into individual genomes and degapped, forming
the sample-wise references. Reads for each run are simulated from each genome
using mason2. This simulator includes a realistic Illumina sequencing error
profile, and we do a couple of "technical replicates" per run. The number of
reads for a run is randomly drawn from a normal distribution with mean of the
average coverage and a coefficient of variation of 0.3. This is roughly what
you'd get in real world sequencing experiments. We then use my trimit tool to
remove read-through and poor quality regions using a windowed trimming
algorithm. This will remove simulated sequences that would not be tolerated in
real datasets due to high error rate or adaptor contamination.

These reads are hashed with khmer into a hash with a tablesize of 5e7. Again
this is to reduce the computational time of the simulation and should probably
be increased for the final version of this simulation. We use k=20, and disable
bigcount (counts past 255) as we have no way of utilising this feature.

We then run kWIP over each set of runs, in both weighted and unweighted mode.
The kernel and distance matrices are plotted using the img.R script from the
kWIP repo.

We also calculate a distance matrix directly from the original tree between
samples. Then, we calculate a least-squares regression and use the sum squared
residuals as the measure of fit.

Possible other measures of fit:

  - NJ or hclust dendrogram then Robertson-Foulds distance?
  - Matrix similarity measure? (Mantel/some other??)


Results
^^^^^^^
