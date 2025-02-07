# base image: Ubuntu
FROM ubuntu:24.10

ENV DEBIAN_FRONTEND=noninteractive

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

# install perl, python3, blast and braker3
RUN apt-get update --fix-missing && \
apt-get install -y wget curl gcc gawk grep psmisc sed coreutils unzip openjdk-11-jdk g++ make perl cpanminus git libboost-all-dev libmysql++-dev libsqlite3-dev libgsl-dev liblpsolve55-dev python3 cmake ncbi-blast+ libbz2-dev libclang-dev liblzma-dev bzip2 libncurses5-dev libncursesw5-dev libglib2.0-dev sudo python3-biopython bamtools libbamtools-dev samtools libhts-dev && \
cpanm Scalar::Util::Numeric File::HomeDir Parallel::ForkManager Math::Utils YAML Hash::Merge MCE::Mutex Cwd Data::Dumper File::Path File::Spec File::Temp FindBin Getopt::Long Hash::Merge List::Util MCE::Mutex Math::Utils Parallel::ForkManager Statistics::LineFit Storable Thread::Queue YAML::XS threads && \
wget https://github.com/Gaius-Augustus/BRAKER/archive/refs/tags/v3.0.8.tar.gz && \
tar -xvf v3.0.8.tar.gz

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

# install GeneMark-ETP
# there is an error with selectSupportedSubsets.py script, which is ModuleNotFoundError: No module named 'statistics.predictionAnalysis'; 'statistics' is not a package
# the developer recommends downloading this script via git (https://github.com/Gaius-Augustus/BRAKER/tree/report/scripts), along with three other scripts, and replacing them.
# Git issue : https://github.com/Gaius-Augustus/BRAKER/issues/358
RUN git clone https://github.com/gatech-genemark/GeneMark-ETP && \
cp /GeneMark-ETP/tools/* /usr/local/bin/ && \
cp -R /GeneMark-ETP/bin/* /usr/local/bin/ && \
cp -R /GeneMark-ETP/bin/gmes/* /usr/local/bin/ && \
cp -R /GeneMark-ETP/bin/gmst/* /usr/local/bin/ && \
cp -R /GeneMark-ETP/bin/GeneMarkSTFiltering/* /usr/local/bin/
COPY predictionAnalysis_for_GeneMark-ETP predictionAnalysis/
RUN cp /predictionAnalysis/* /usr/local/bin/

# Install Augustus
RUN git clone https://github.com/Gaius-Augustus/Augustus.git && \
cd /Augustus && \
make augustus && \
cd /Augustus/auxprogs/joingenes/ && \
make && \
cd /Augustus/auxprogs/bam2hints/ && \
make && \
cd /Augustus/auxprogs/aln2wig/ && \
make && \
cd /Augustus/auxprogs/bam2wig/ && \
make && \
cd /Augustus/auxprogs/compileSpliceCands/ && \
make && \
cd /Augustus/auxprogs/filterBam/ && \
make && \
cd /Augustus/auxprogs/homGeneMapping/ && \
make && \
cd /Augustus/auxprogs/utrrnaseq/ && \
make && \
ln -s /Augustus/auxprogs/utrrnaseq/utrrnaseq /Augustus/bin/ && \
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

# install TSEBRA
RUN git clone https://github.com/Gaius-Augustus/TSEBRA && \
cp /TSEBRA/bin/*  /usr/local/bin/
