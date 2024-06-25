# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

# install perl and diamond
RUN apt-get update --fix-missing && \
apt-get install -y wget gcc g++ cmake libdb-dev make perl cpanminus r-base libbio-perl-perl sudo && \
cpanm --force Statistics::R Clone Graph::Directed LWP::UserAgent Carp Sort::Naturally File::Share File::ShareDir::Install Moose YAML LWP::Protocol::https Term::ProgressBar
RUN wget https://github.com/bbuchfink/diamond/archive/refs/tags/v2.0.15.tar.gz && \
tar -xvf v2.0.15.tar.gz && \
cd diamond-2.0.15 && \
mkdir bin && \
cd bin && \
cmake .. && \
make -j4 && \
sudo make install

# install AGAT
RUN wget https://github.com/NBISweden/AGAT/archive/refs/tags/v1.4.0.tar.gz && \
tar -xvf v1.4.0.tar.gz && \
cd AGAT-1.4.0 && \
perl Makefile.PL && \
make && \
make install