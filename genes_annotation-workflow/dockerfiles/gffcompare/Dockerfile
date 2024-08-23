# base image: Ubuntu
FROM ubuntu:24.10

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

# install gffcompare
RUN apt-get update --fix-missing && \
apt-get install -y wget gcc g++ make && \
wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.12.6.tar.gz && \
tar -xvf gffcompare-0.12.6.tar.gz && \
cd gffcompare-0.12.6/ && \
make && \
ln -s gffcompare /usr/local/bin/ && \
ln -s trmap /usr/local/bin/
