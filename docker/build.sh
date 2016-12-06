#!/bin/bash

set -xeuo pipefail

# Tool versions
declare -A VERS
VERS[snakemake]=3.8.2
VERS[khmer]=2.0
VERS[skbio]=0.5.1
VERS[ete3]=3.0.0b35
VERS[scrm]=1.7.2
VERS[mason2]=2.0.5
VERS[dawg]=f9ebcd8cff   # This commit was the dev HEAD for initial experiments
VERS[trimit]=0.2.5
VERS[kwip]=0.2.0
VERS[mash]=1.1.1

## More about DAWG version:
# We require features from the develop branch of DAWG, which will become dawg2.
# While it was not when these experiments were written, this tool is now under
# active development, making it no longer possible to simply use the develop
# branch reproducibly. Therefore, we use the state that the develop branch was
# in when these experiments started.


################################################################################
#                                 Apt packages                                 #
################################################################################

sed -i -e 's/archive.ubuntu.com/au.archive.ubuntu.com/' /etc/apt/sources.list
apt-get update -yy
apt-get upgrade -yy

apt-get install -yy              \
    vim                          \
    curl                         \
    build-essential              \
    cmake                        \
    pkg-config                   \
    libboost-dev                 \
    libboost-program-options-dev \
    libbz2-dev                   \
    libgsl-dev                   \
    zlib1g-dev                   \
    cython3                      \
    python3-dev                  \
    python3-pip                  \
    python3-numpy                \
    python3-scipy                \
    python3-docopt               \
    python3-matplotlib           \
    python3-yaml                 \
    python3-six                  \
    python3-tz                   \
    python3-dateutil             \



################################################################################
#                                     Pip                                      #
################################################################################

pip3 install --pre                \
    snakemake==${VERS[snakemake]} \
    khmer==${VERS[khmer]}         \
    scikit-bio==${VERS[skbio]}    \
    ete3==${VERS[ete3]}           \


################################################################################
#                            Source/binary tarballs                            #
################################################################################

cd /usr/local/src

##########
#  SCRM  #
##########

# Downloads static binary directly
curl -LS -o /usr/local/bin/scrm \
    https://github.com/scrm/scrm/releases/download/v${VERS[scrm]}/scrm-x64-static
chmod +x /usr/local/bin/scrm


############
#  Mason2  #
############

curl -LSO http://packages.seqan.de/mason2/mason2-${VERS[mason2]}-Linux-x86_64.tar.xz
tar xvf mason2-${VERS[mason2]}-Linux-x86_64.tar.xz
mv mason2-${VERS[mason2]}-Linux-x86_64/bin/* /usr/local/bin
rm -rf /usr/local/src/*


##############################
#  Dawg Development Release  #
##############################

tarname=dawg_${VERS[dawg]}.tar.gz
# Install dawg2
curl -LS -o ${tarname}  \
    https://github.com/reedacartwright/dawg/archive/${VERS[dawg]}.tar.gz
tar xvf ${tarname}
cd dawg-${VERS[dawg]}*/
mkdir -p build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=..
make all install
mv ../bin/dawg /usr/local/bin/dawg2
mv ../lib/* /usr/local/lib
cd ../..
rm -rf /usr/local/src/*
unset tarname

############
#  Trimit  #
############

curl -LSO \
    https://github.com/kdmurray91/libqcpp/releases/download/${VERS[trimit]}/trimit_${VERS[trimit]}_amd64.tar.gz

# extracts bin/trimit to /usr/local
tar xvf trimit_${VERS[trimit]}_amd64.tar.gz -C /usr/local --strip-components=1


##########
#  kWIP  #
##########

tarname=kWIP_${VERS[kwip]}.tar.gz
curl -LS -o $tarname \
    https://github.com/kdmurray91/kWIP/archive/${VERS[kwip]}.tar.gz
tar xvf $tarname
cd kWIP-${VERS[kwip]}/
mkdir build && cd build
cmake ..
make all install
make test
rm -rf /usr/local/src/*

##########
#  MASH  #
##########

curl -LS \
    https://github.com/marbl/Mash/releases/download/v${VERS[mash]}/mash-Linux64-v${VERS[mash]}.tar.gz \
    | tar xz --strip-components 1 -C /usr/local/bin
rm -rf /usr/local/src/*


################################################################################
#                                   Cleanup                                    #
################################################################################

apt-get purge -yy   \
    curl            \
    build-essential \
    cmake           \
    pkg-config      \
    libboost-dev    \
    libbz2-dev      \
    zlib1g-dev      \
    python3-dev     \

apt-get autoremove -y
apt-get clean -y
rm -rf /var/cache/apt/archives/* /var/lib/apt/lists/* /usr/share/locale/*
rm -rf /root/.cache
