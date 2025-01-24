# Use Ubuntu 20.04 LTS as the base image
FROM ubuntu:20.04

# Set environment variables to prevent Python from writing .pyc files to disk
ENV PYTHONUNBUFFERED=1

# Install necessary system packages
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    build-essential \
    libopenblas-dev \
    liblapack-dev \
    libatlas-base-dev \
    gfortran \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip and setuptools
RUN pip3 install --no-cache-dir --upgrade pip setuptools wheel

# Install the specified Python packages
RUN pip3 install --no-cache-dir numpy==1.21.0

RUN pip3 install --no-cache-dir \
    kaleido==0.2.1 \
    pandas==1.1.5 \
    plotly==5.18.0 \
    biopython==1.79 \
    matplotlib==3.5.2 \
    scipy==1.10.1 \
    colorlover==0.3.0 \
    tqdm==4.62.3 \
    networkx==3.0b1

# Copy the aegis directory into the image
COPY aegis /aegis

# Set the working directory
WORKDIR /aegis

# Add the aegis directory to PYTHONPATH
ENV PYTHONPATH="/aegis"

# Default command to keep the container running in interactive mode
CMD ["/bin/bash"]

