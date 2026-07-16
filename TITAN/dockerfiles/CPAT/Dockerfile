# base image: Rocky linux
FROM rockylinux@sha256:d644d203142cd5b54ad2a83a203e1dee68af2229f8fe32f52a30c6e1d3c3a9e0

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

# Install base utilities
RUN dnf -y upgrade && \
dnf install -y wget java-11-openjdk unzip git bash gawk grep psmisc sed gcc gcc-c++ make zlib-devel python3 zlib-devel ncurses-devel bzip2 xz-devel glib2-devel python3-numpy python3-pip dnf-plugins-core && \
dnf -y config-manager --set-enabled crb && \
dnf -y install epel-release && \
dnf -y install R

# Install CPAT
RUN pip3 install pysam 
RUN pip3 install CPAT

# command : cpat -h
# make_hexamer_tab -h
# make_logitModel -h
