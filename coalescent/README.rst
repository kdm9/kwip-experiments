Coalescent Simulation
=====================


To re-run this simulation in the docker container: ::

    # Build the image
    docker build --tag kwip-expt:latest \
        https://github.com/kdmurray91/kwip-experiments-docker.git
    # Check that builds fine, you should see "Successfully build image" or
    # similar.

    # Create a new container (aka VM) to run simulations within
    docker run -it --name kwip-container kwip-expt
    # you should now be inside a container with all required software

    # Get the simulation repo
    git clone https://github.com/kdmurray91/2015_kwip-coalescent coalescent
    cd coalescent

    # Run the simulation with NCPU cpus. Change NCPU to the number of jobs you
    # want concurrently run
    snakemake -j <NCPU>


To work with docker on OSX, see `this guide <https://docs.docker.com/mac/>`_.

