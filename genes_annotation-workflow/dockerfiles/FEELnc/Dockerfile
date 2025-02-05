# base image: Rocky linux
FROM rockylinux:9

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

# Install base utilities
RUN dnf -y upgrade && \
dnf install -y wget java-11-openjdk unzip git bash gawk grep psmisc sed gcc gcc-c++ make git-all python3 zlib-devel ncurses-devel bzip2 xz-devel glib2-devel python which perl cpanminus expat-devel procps dnf-plugins-core && \
cpanm --force Parallel::ForkManager Bio::DB::Fasta Bio::DB::SeqFeature::Store && \
dnf -y config-manager --set-enabled crb && \
dnf -y install epel-release && \
dnf -y install R

# Install FEELnc from Git
RUN git clone https://github.com/tderrien/FEELnc.git
ENV FEELNCPATH="/FEELnc/"
ENV PERL5LIB="$PERL5LIB:/FEELnc/lib/"
ENV PATH="$PATH:/FEELnc/scripts/"
ENV PATH="$PATH:/FEELnc/utils/"
ENV PATH="$PATH:/FEELnc/bin/LINUX/"