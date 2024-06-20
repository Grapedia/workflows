# base image: Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

# install perl and braker2
RUN apt-get update --fix-missing && \
apt-get install -y wget gcc g++ make perl cpanminus && \
cpanm Scalar::Util::Numeric File::HomeDir Parallel::ForkManager Math::Utils YAML && \
wget https://github.com/Gaius-Augustus/BRAKER/archive/refs/tags/v2.1.6.tar.gz && \
tar -xvf v2.1.6.tar.gz

# install spaln
RUN wget https://github.com/ogotoh/spaln/archive/refs/tags/Ver.2.3.3.tar.gz && \
tar -xvf Ver.2.3.3.tar.gz && \
cd spaln-Ver.2.3.3/src/ && \
./configure && \
make && \
make install && \
mv spaln sortgrcd makmdm makdbs /usr/local/bin

# install ProtHint
RUN wget https://github.com/gatech-genemark/ProtHint/releases/download/v2.6.0/ProtHint-2.6.0.tar.gz && \
tar -xvf ProtHint-2.6.0.tar.gz

# install GeneMark
# Link to download : http://topaz.gatech.edu/GeneMark/license_download.cgi
# the license is valid 200 days
# download GeneMark-ES/ET/EP+ manually and copy the file to the dockerfile
# tar -xvf gmes_linux_64_4.tar.gz
# gunzip gm_key_64.gz
COPY gmes_linux_64_4 GeneMark/
COPY gm_key_64 /root/.gm_key
