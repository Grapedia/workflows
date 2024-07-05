# base image: Ubuntu
FROM ubuntu:24.10

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

# install perl, python3, blast and braker2
RUN apt-get update --fix-missing && \
apt-get install -y wget gcc g++ make perl cpanminus git libboost-all-dev libmysql++-dev libsqlite3-dev libgsl-dev liblpsolve55-dev python3 cmake ncbi-blast+ libbz2-dev libclang-dev liblzma-dev bzip2 libncurses5-dev libncursesw5-dev libglib2.0-dev sudo python3-biopython && \
cpanm Scalar::Util::Numeric File::HomeDir Parallel::ForkManager Math::Utils YAML Hash::Merge MCE::Mutex && \
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
# download GeneMark-ES/ET/EP+ manually and copy the file to the dockerfile folder
# tar -xvf gmes_linux_64_4.tar.gz
# gunzip gm_key_64.gz
COPY gmes_linux_64_4 GeneMark/
COPY gm_key_64 /root/.gm_key

# Install Augustus
RUN git clone https://github.com/Gaius-Augustus/Augustus.git && \
cd /Augustus && \
make augustus && \
cd /Augustus/auxprogs/joingenes/ && \
make && \
chmod -R 777 /Augustus

# Install bamtools
RUN wget https://github.com/pezmaster31/bamtools/archive/refs/tags/v2.5.2.tar.gz && \
tar -xvf v2.5.2.tar.gz && \
cd bamtools-2.5.2/ && \
mkdir build && \
cd build && \
cmake -DCMAKE_INSTALL_PREFIX=/usr/local/bin/ .. && \
make && \
make install && \
ln -s /usr/local/bin/bin/bamtools /usr/local/bin/bamtools

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
tar -xjf samtools-1.9.tar.bz2 && \
rm samtools-1.9.tar.bz2 && \
cd samtools-1.9 && \
./configure && \
make && \
make install

# install diamond
RUN wget https://github.com/bbuchfink/diamond/archive/refs/tags/v2.1.9.tar.gz && \
tar -xvf v2.1.9.tar.gz && \
cd diamond-2.1.9 && \
mkdir bin && \
cd bin && \
cmake .. && \
make -j4 && \
sudo make install

# install cdbfasta
RUN git clone https://github.com/gpertea/cdbfasta.git && \
cd cdbfasta/ && \
make && \
cp cdbfasta cdbyank /usr/local/bin/
