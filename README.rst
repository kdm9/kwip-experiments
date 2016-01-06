==============================
kWIP experiments Docker image
==============================

Usage
^^^^^

::
    docker build --tag kwip-expt:latest https://github.com/kdmurray91/kwip-experiments-docker.git
    docker run -it --name kwip-container kwip-expt


Software
^^^^^^^^


All software needed to run the experiments should be included, if something isn't, raise an issue.

This contains the following software:

 - kWIP
 - libqcpp/trimit
 - dawg 2.0+
 - scrm
 - seqan's mason2 simulator
 - scientific python stack (np, mpl, scipy, pandas, jupyter etc)
 - scikit-bio
 - pycogent
 - ete
 - snakemake
 
