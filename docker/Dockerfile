FROM debian:jessie
MAINTAINER Kevin Murray <spam@kdmurray.id.au>

# Install debian packages
RUN sed -i -e 's/httpredir.debian.org/mirrors.kernel.org/' /etc/apt/sources.list \
    && apt-get update \
    && apt-get upgrade -yy \
    && apt-get install -yy \
        less                    \
        ncdu                    \
        screen                  \
        tree                    \
        vim                     \
        wget                    \
        bison                   \
        build-essential         \
        gfortran                \
        cmake                   \
        flex                    \
        git                     \
        pkg-config	        \
        libboost-dev            \
        libboost-program-options-dev \
        libbz2-dev              \
        libgsl0-dev             \
        libyaml-cpp-dev         \
        libxml2-dev             \
        libxslt1-dev            \
        libblas-dev             \
        libfreetype6-dev        \
        libpng12-dev            \
        liblapack-dev           \
        libsuitesparse-dev      \
        zlib1g-dev              \
        cython3                 \
        python3-dev             \
        python3-pip             \
    && apt-get autoremove       \
    && apt-get clean -y         \
    && rm -rf /var/cache/apt/archives/* /var/lib/apt/lists/*

# Install python applications & pacakges
RUN pip3 install cython \
    && pip3 install numpy>=1.9.2 \
    && pip3 install scipy
RUN pip3 install --pre ete3 \
    && pip3 install snakemake khmer==2.0 scikit-bio \
                    matplotlib docopt

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
ADD https://github.com/kdmurray91/kWIP/releases/download/0.2.0-rc2/kwip-binaries_0.2.0-rc2.tar.xz /usr/local/src/
RUN cd /usr/local/ && \
    tar xvf src/kwip-binaries_*.tar.xz && \
    rm -rf src/*.tar*

# Obtain the experiments
WORKDIR /root
RUN git clone https://github.com/kdmurray91/kwip-experiments
