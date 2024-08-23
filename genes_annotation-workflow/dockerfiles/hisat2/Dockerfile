# base image: Ubuntu
FROM ubuntu:24.10

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

# install hisat2-2.2.1
RUN apt-get update --fix-missing && \
apt-get install -y wget gcc g++ make unzip python3 samtools && \
ln -s /usr/bin/python3 /usr/bin/python  && \
wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download && \
unzip download