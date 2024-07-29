# base image: Rocky linux
FROM rockylinux:9

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

# Install base utilities
RUN dnf -y upgrade && \
dnf install -y wget gcc gcc-c++ make git-all zlib-devel python3 zlib-devel ncurses-devel bzip2 bzip2-devel xz-devel glib2-devel python which

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
tar -xjf samtools-1.9.tar.bz2 && \
cd samtools-1.9 && \
./configure && \
make && \
make install

# install bedtools
RUN wget -O bedtools.tar.gz https://github.com/arq5x/bedtools2/archive/refs/tags/v2.30.0.tar.gz && \
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
