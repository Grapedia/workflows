# base image: Ubuntu
FROM ubuntu:24.10

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

# install perl, python3, blast and braker3
RUN apt-get update --fix-missing && \
apt-get install -y curl wget gcc g++ make perl cpanminus git libboost-all-dev libmysql++-dev libsqlite3-dev libgsl-dev liblpsolve55-dev python3 cmake ncbi-blast+ libbz2-dev libclang-dev liblzma-dev bzip2 libncurses5-dev libncursesw5-dev libglib2.0-dev sudo python3-biopython bamtools libbamtools-dev libhts-dev && \
curl -L https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 | tar -jxvf - && \
ln -s /minimap2-2.28_x64-linux/minimap2 /usr/local/bin/minimap2

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
tar -xjf samtools-1.9.tar.bz2 && \
rm samtools-1.9.tar.bz2 && \
cd samtools-1.9 && \
./configure && \
make && \
make install