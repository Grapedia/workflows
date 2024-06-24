# base image: debian
FROM debian:unstable-slim AS builder

# File Author / Maintainer
MAINTAINER Amandine Velt <amandine.velt@inrae.fr>

RUN apt-get update && \
    apt-get install --no-install-recommends -y \
        curl libc6 \
        ca-certificates \
        build-essential

RUN curl -L https://github.com/guigolab/geneid/archive/v1.4.4.tar.gz | \
    tar xz && \
    cd geneid-1.4.4 && \
    make BIN=/build && \
    mv /build/geneid /usr/local/bin/
