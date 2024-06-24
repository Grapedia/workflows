# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

# Install base utilities
RUN apt-get update && \
apt-get install -y build-essential && \
apt-get install -y wget gcc g++ make

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
/bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

RUN eval "$(conda shell.bash hook)" && \
conda create -n EDTA && \
conda activate EDTA && \
conda install -y -n base --override-channels -c conda-forge mamba 'python_abi=*=*cp*' && \
mamba install --yes -c conda-forge -c bioconda edta