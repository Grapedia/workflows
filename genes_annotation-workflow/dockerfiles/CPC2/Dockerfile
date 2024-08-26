# base image: Ubuntu
FROM ubuntu:24.10

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

# install CPC2
RUN apt-get update --fix-missing && \
apt-get install -y wget gcc g++ make python3 
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN apt-get install -y python3-biopython python3-six
RUN wget https://github.com/gao-lab/CPC2_standalone/archive/refs/tags/v1.0.1.tar.gz && \
tar -xvf v1.0.1.tar.gz && \
cd CPC2_standalone-1.0.1/ && \
export CPC_HOME="/CPC2_standalone-1.0.1/" && \
cd libs/libsvm && \
gzip -dc libsvm-3.18.tar.gz | tar xf - && \
cd libsvm-3.18 && \
make clean && \
make

# command line example : cd /CPC2_standalone-1.0.1/; ./bin/CPC2.py -h
