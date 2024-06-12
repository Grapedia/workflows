# base image: Ubuntu
FROM python:3

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

# install gffutils
RUN apt-get update --fix-missing && \
apt-get install -y wget gcc g++ make zlib1g zlib1g-dev software-properties-common
RUN pip install gffutils

# install GlimmerHMM
RUN wget https://ccb.jhu.edu/software/glimmerhmm/dl/GlimmerHMM-3.0.4.tar.gz && \
tar -xzvf GlimmerHMM-3.0.4.tar.gz && \
cd GlimmerHMM && \
mv bin/* /usr/local/bin && \
mv train/* /usr/local/bin
