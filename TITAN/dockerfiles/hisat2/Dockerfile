# base image: Ubuntu
FROM ubuntu@sha256:f995e05e8adc3292853cc37e6edda72351f8002ce7469a29322d19e01529cb9f

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

# install hisat2-2.2.1
RUN apt-get update --fix-missing && \
apt-get install -y wget gcc g++ make unzip python3 samtools curl gawk grep bash psmisc sed coreutils unzip openjdk-11-jdk git && \
ln -s /usr/bin/python3 /usr/bin/python  && \
wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download && \
unzip download