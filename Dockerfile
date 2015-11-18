FROM debian:testing
MAINTAINER Kevin Murray <spam@kdmurray.id.au>

RUN apt-get update && apt-get upgrade -yy
RUN apt-get install -yy \
        python-pip      \
        python3-pip     \
        python-lxml     \
        python-matplotlib \
        python-numpy    \
        python-pandas   \
        python-scipy    \
        python-six      \
        cython


# Install python applications & pacakges
RUN pip3 install snakemake khmer==2.0
RUN pip install ete2 scikit-bio

# Install scrm
ADD https://github.com/scrm/scrm/releases/download/v1.6.1/scrm-x64-static.tar.gz /usr/local/src
RUN cd /usr/local/src && \
    tar xvf scrm-x64-static.tar.gz && \
    mv scrm /usr/local/bin && \
    rm -rf doc/manual.html scrm-x64-static.tar.gz

# Install mason2
ADD http://packages.seqan.de/mason2/mason2-2.0.1-Linux-x86_64.tar.bz2 /usr/local/src
RUN cd /usr/local/src && \
    tar xvf mason2-2.0.1-Linux-x86_64.tar.bz2 && \
    mv mason2-2.0.1-Linux-x86_64/bin/* /usr/local/bin && \
    rm -rf mason2-2.0.1-Linux-x86_64 mason2-2.0.1-Linux-x86_64.tar.bz2

