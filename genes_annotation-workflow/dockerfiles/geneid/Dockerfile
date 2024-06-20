# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

# install geneid
RUN apt-get update --fix-missing && \
apt-get install -y wget gcc g++ make zlib1g zlib1g-dev python libbz2-dev libclang-dev liblzma-dev bzip2 libncurses5-dev libncursesw5-dev libglib2.0-dev && \
wget -O bedtools.tar.gz https://github.com/arq5x/bedtools2/archive/refs/tags/v2.30.0.tar.gz && \
tar -xzvf bedtools.tar.gz && \
cd bedtools2-2.30.0 && \
make && \
mv bin/* /usr/local/bin

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
tar -xjf samtools-1.9.tar.bz2 && \
rm samtools-1.9.tar.bz2 && \
cd samtools-1.9 && \
./configure && \
make && \
make install

# install exonerate
RUN wget https://github.com/nathanweeks/exonerate/archive/refs/tags/v2.4.0.tar.gz  && \
tar -xvf v2.4.0.tar.gz && \
cd exonerate-2.4.0 && \
./configure && \
make && \
make install
