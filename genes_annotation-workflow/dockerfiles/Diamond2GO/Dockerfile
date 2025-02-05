# base image: Rocky linux
FROM rockylinux:9

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

# Install base utilities
RUN dnf -y upgrade && \
dnf install -y --nobest wget unzip sed psmisc grep gawk java-11-openjdk gcc bash gcc-c++ make git git-all zlib-devel python3 zlib-devel ncurses-devel bzip2 xz-devel glib2-devel python which

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
/bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

# install dependencies : GIT-LFS, perl, bioperl, Diamond, Interproscan and some perl modules
RUN eval "$(conda shell.bash hook)" && \
conda config --add channels conda-forge && \
conda config --set channel_priority strict && \
conda install git-lfs && \
conda install bioconda::perl-bioperl && \
conda install conda-forge::perl-app-cpanminus && \
cpanm install Getopt::Std && \
conda install gxx_linux-64 && \
cpanm install Scalar::Util && \
conda install bioconda::diamond && \
conda install bioconda::interproscan && \
cpanm install List::AllUtils && \
cpanm install LWP && \
cpanm install LWP::Protocol::https && \
cpanm install Mozilla::CA && \
cpanm install Time::HiRes && \
cpanm install XML::Simple

# download of Diamond2GO
RUN git-lfs clone https://github.com/rhysf/Diamond2GO.git && \
cd /Diamond2GO && \
git lfs install && \
git-lfs pull && \
chmod 755 Diamond2go.pl

# example of command without interproscan
# the default NCBI non-redundant database was downloaded on the 14th May 2023
# perl /Diamond2GO/Diamond2go.pl -d /Diamond2GO/resources/nr_clean_d2go.dmnd \
# -q proteins.fasta -t protein

# example of command with interproscan only on genes/proteins with no diamond2go hits
# perl /Diamond2GO/Diamond2go.pl -d /Diamond2GO/resources/nr_clean_d2go.dmnd \
# -q proteins.fasta -t protein -s 1234 -i h -z email@email.com
