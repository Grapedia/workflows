# base image: micromamba 1.5.10
FROM mambaorg/micromamba@sha256:e3797091302382ea841498bc93a7b0a50f7c1448333d5e946d2d1608d0c5f43d

# File Author / Maintainer
LABEL maintainer="Amandine Velt <amandine.velt@inrae.fr>"

# Mender (split-gene detection and merging, https://github.com/simrbase/mender) is a set of Perl
# scripts with no release or package: pinned here to a commit. Its dependencies come from bioconda.
ARG MENDER_COMMIT=959fabbbc93f0e14180488262eea74b9394dd6f4

USER root
RUN micromamba install -y -n base -c conda-forge -c bioconda \
      perl=5.32 diamond=2.1.10 gffread=0.12.7 genometools-genometools=1.6.5 bedtools=2.31.1 \
      mafft=7.526 samtools=1.21 python=3.12 git procps-ng && \
    /opt/conda/bin/git clone https://github.com/simrbase/mender.git /opt/mender && \
    /opt/conda/bin/git -C /opt/mender checkout ${MENDER_COMMIT} && \
    echo ${MENDER_COMMIT} > /opt/mender/COMMIT && \
    rm -rf /opt/mender/.git /opt/mender/*.png && \
    chmod +x /opt/mender/*.pl /opt/mender/*.sh && \
    micromamba clean -a -y

ENV PATH=/opt/conda/bin:/opt/mender:$PATH
ENV MENDER_DIR=/opt/mender
