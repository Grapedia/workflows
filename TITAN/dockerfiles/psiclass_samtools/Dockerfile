# base image: Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

ENV DEBIAN_FRONTEND=noninteractive

# install PsiCLASS-1.0.2
RUN apt-get update --fix-missing && \
apt-get install -y curl gawk grep psmisc sed coreutils unzip openjdk-11-jdk git wget gcc g++ make zlib1g zlib1g-dev python libbz2-dev libclang-dev liblzma-dev bzip2 libncurses5-dev libncursesw5-dev libglib2.0-dev && \
wget https://github.com/splicebox/PsiCLASS/archive/refs/tags/v1.0.2.tar.gz && \
tar -xvf v1.0.2.tar.gz && \
cd PsiCLASS-1.0.2 && \
make
