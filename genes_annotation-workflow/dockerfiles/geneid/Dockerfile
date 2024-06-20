# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

# install geneid
RUN apt-get update --fix-missing && \
apt-get install -y wget gcc g++ make zlib1g zlib1g-dev python libbz2-dev libclang-dev liblzma-dev bzip2 libncurses5-dev libncursesw5-dev libglib2.0-dev && \
wget https://genome.crg.es/pub/software/geneid/geneid_v1.4.4.Jan_13_2011.tar.gz && \
tar -xzvf geneid_v1.4.4.Jan_13_2011.tar.gz && \
cd geneid && \
make && \
mv bin/* /usr/local/bin
