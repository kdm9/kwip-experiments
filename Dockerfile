FROM debian:testing
MAINTAINER Kevin Murray <spam@kdmurray.id.au>

# Install debian packages
RUN apt-get update && apt-get upgrade -yy

## Base applications
RUN apt-get install -yy --no-install-recommends \
        htop            \
        less            \
        ncdu            \
        screen          \
        tree            \
        vim             \
        wget

## Base devl stuff
RUN apt-get install -yy --no-install-recommends \
        bison           \
        build-essential \
        cmake           \
        cython          \
        flex            \
        git             \
        libboost-all-dev \
        libbz2-dev      \
        libgsl-dev      \
        libyaml-cpp-dev \
        python-dev      \
        zlib1g-dev


# Python software
RUN apt-get install -yy --no-install-recommends \
        python-cogent   \
        python-docopt   \
        python-lxml     \
        python-matplotlib \
        python-numpy    \
        python-pandas   \
        python-pip      \
        python-scipy    \
        python-six      \
        python-skbio    \
        snakemake

# this is an optional dep of ete2, and it makes the image massive
#        python-qt4      \

# Install python applications & pacakges
RUN pip install ete2 khmer==2.0

# Install scrm
RUN cd /usr/local/src && \
    wget 'https://github.com/scrm/scrm/releases/download/v1.6.1/scrm-x64-static.tar.gz' && \
    tar xvf scrm-x64-static.tar.gz && \
    mv scrm /usr/local/bin && \
    rm -rf doc scrm-x64-static.tar.gz

# Install mason2
RUN cd /usr/local/src && \
    wget 'http://packages.seqan.de/mason2/mason2-2.0.1-Linux-x86_64.tar.bz2' && \
    tar xvf mason2-2.0.1-Linux-x86_64.tar.bz2 && \
    mv mason2-2.0.1-Linux-x86_64/bin/* /usr/local/bin && \
    rm -rf mason2-2.0.1-Linux-x86_64 mason2-2.0.1-Linux-x86_64.tar.bz2

# Install dawg2
RUN cd /usr/local/src && \
    git clone https://github.com/reedacartwright/dawg/ dawg && \
    cd dawg && \
    git checkout develop && \
    mkdir -p build && cd build && \
    cmake .. && \
    make && \
    make install && \
    mv /usr/local/bin/dawg /usr/local/bin/dawg2 && \
    cd ../.. && \
    rm -rf dawg

# Install libqc++ for trimit
RUN cd /usr/local/src && \
    git clone https://github.com/kdmurray91/libqcpp && \
    cd libqcpp && \
    mkdir -p build && cd build && \
    cmake .. && \
    make && \
    make install && \
    cd ../.. && \
    rm -rf libqcpp

# Install kWIP
RUN cd /usr/local/src && \
    git clone https://github.com/kdmurray91/kwip && \
    cd kwip && \
    mkdir -p build && cd build && \
    cmake .. && \
    make all install && \
    cd ../.. && \
    rm -rf kwip
