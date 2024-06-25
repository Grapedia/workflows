# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

# install perl and gffread
RUN apt-get update --fix-missing && \
apt-get install -y wget gcc g++ make zlib1g zlib1g-dev perl cpanminus && \
cpanm URI DBI && \
wget https://github.com/gpertea/gffread/releases/download/v0.12.7/gffread-0.12.7.Linux_x86_64.tar.gz && \
tar -xvf gffread-0.12.7.Linux_x86_64.tar.gz && \
cd gffread-0.12.7.Linux_x86_64/ && \
mv gffread /usr/local/bin

# install evidencemodeler
RUN wget https://github.com/EVidenceModeler/EVidenceModeler/archive/refs/tags/v1.1.1.tar.gz && \
tar -xvf v1.1.1.tar.gz && \
cd EVidenceModeler-1.1.1 && \
mv evidence_modeler.pl /usr/local/bin/ && \
mv EvmUtils/*.pl /usr/local/bin/

ENV PERL5LIB="${PERL5LIB}:/EVidenceModeler-1.1.1/PerlLib/"