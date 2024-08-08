# base image: Ubuntu
FROM ubuntu:24.10

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

# install perl, python3, blast and braker2
RUN apt-get update --fix-missing && \
apt-get install -y wget gcc g++ make libboost-all-dev libmysql++-dev libsqlite3-dev libgsl-dev liblpsolve55-dev python3 cmake ncbi-blast+ libbz2-dev libclang-dev liblzma-dev bzip2 libncurses5-dev libncursesw5-dev libglib2.0-dev sudo

# install stringtie
RUN wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.2.3.Linux_x86_64.tar.gz && \
tar -xvf stringtie-2.2.3.Linux_x86_64.tar.gz && \
cd stringtie-2.2.3.Linux_x86_64 && \
mv stringtie /usr/local/bin
