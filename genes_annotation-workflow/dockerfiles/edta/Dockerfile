# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

# Install base utilities
RUN apt-get update && \
apt-get install -y build-essential && \
apt-get install -y wget gcc g++ make git zlib1g zlib1g-dev python libbz2-dev libclang-dev liblzma-dev bzip2 libncurses5-dev libncursesw5-dev libglib2.0-dev

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
tar -xjf samtools-1.9.tar.bz2 && \
rm samtools-1.9.tar.bz2 && \
cd samtools-1.9 && \
./configure && \
make && \
make install

# install bedtools
RUN apt-get update --fix-missing && \
apt-get install -y wget gcc g++ make zlib1g zlib1g-dev python libbz2-dev libclang-dev liblzma-dev bzip2 libncurses5-dev libncursesw5-dev libglib2.0-dev && \
wget -O bedtools.tar.gz https://github.com/arq5x/bedtools2/archive/refs/tags/v2.30.0.tar.gz && \
tar -xzvf bedtools.tar.gz && \
cd bedtools2-2.30.0 && \
make && \
mv bin/* /usr/local/bin

# Retrieve EDTA from GIT
RUN git clone https://github.com/oushujun/EDTA.git

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
/bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

RUN eval "$(conda shell.bash hook)" && \
cd /EDTA  && \
conda install -y -n base --override-channels -c conda-forge mamba 'python_abi=*=*cp*' && \
mamba env create -f EDTA_2.2.x.yml && \
conda activate EDTA2
