# base image: Rocky linux
FROM rockylinux:9

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

# Install base utilities
RUN dnf -y upgrade && \
dnf install -y wget gcc gcc-c++ make zlib-devel python3 zlib-devel ncurses-devel bzip2 bzip2-devel xz-devel glib2-devel python3-numpy python3-pip dnf-plugins-core && \
dnf -y config-manager --set-enabled crb && \
dnf -y install epel-release && \
dnf -y install R

# Install CPAT
RUN pip3 install pysam 
RUN pip3 install CPAT

# command : cpat -h
# make_hexamer_tab -h
# make_logitModel -h