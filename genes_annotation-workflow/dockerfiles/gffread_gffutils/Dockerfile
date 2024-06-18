# base image: Ubuntu
FROM python:3

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

# install gffutils
RUN apt-get update --fix-missing && \
apt-get install -y wget gcc g++ make zlib1g zlib1g-dev software-properties-common
RUN pip install gffutils

# install gffread 0.12.7
RUN wget https://github.com/gpertea/gffread/releases/download/v0.12.7/gffread-0.12.7.Linux_x86_64.tar.gz && \
tar -xvf gffread-0.12.7.Linux_x86_64.tar.gz && \
cd gffread-0.12.7.Linux_x86_64/ && \
mv gffread /usr/local/bin
