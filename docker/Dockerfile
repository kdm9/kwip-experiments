FROM debian:testing
MAINTAINER Kevin Murray <spam@kdmurray.id.au>

# Install debian packages
RUN sed -i -e 's/httpredir.debian.org/mirrors.kernel.org/' /etc/apt/sources.list
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

RUN apt-get clean -y

# this is an optional dep of ete2, and it makes the image massive
#        python-qt4      \

# Install python applications & pacakges
RUN pip install ete2 khmer==2.0

# Install scrm
ADD https://github.com/scrm/scrm/releases/download/v1.6.1/scrm-x64-static.tar.gz /usr/local/src/
RUN cd /usr/local/src && \
    tar xvf scrm-x64-static.tar.gz && \
    mv scrm /usr/local/bin && \
    rm -rf doc *.tar*

# Install mason2
ADD http://packages.seqan.de/mason2/mason2-2.0.1-Linux-x86_64.tar.bz2 /usr/local/src/
RUN cd /usr/local/src && \
    tar xvf mason2-*-Linux-x86_64.tar.* && \
    mv mason2-*-Linux-x86_64/bin/* /usr/local/bin && \
    rm -rf mason2-* *.tar*

# Install dawg2
ADD https://github.com/reedacartwright/dawg/archive/develop.tar.gz /usr/local/src/
RUN cd /usr/local/src && \
    tar xvf develop.tar.gz && \
    cd dawg-develop && \
    mkdir -p build && cd build && \
    cmake .. && \
    make && \
    make install && \
    mv /usr/local/bin/dawg /usr/local/bin/dawg2 && \
    cd ../.. && \
    rm -rf dawg* *.tar*

# Install libqc++ for trimit
ADD https://github.com/kdmurray91/libqcpp/archive/0.1.9.tar.gz /usr/local/src/
RUN cd /usr/local/src && \
    tar xvf 0.1.9.tar.gz && \
    cd libqcpp-* && \
    mkdir -p build && cd build && \
    cmake .. && \
    make && \
    make install && \
    cd ../.. && \
    rm -rf libqcpp* *.tar*

# Install kWIP
RUN cd /usr/local/src && \
    git clone https://github.com/kdmurray91/kwip && \
    cd kwip && \
    mkdir -p build && cd build && \
    cmake .. && \
    make all install && \
    cd ../.. && \
    rm -rf kwip
