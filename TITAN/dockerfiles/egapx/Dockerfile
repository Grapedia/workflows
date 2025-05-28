FROM ncbi/egapx:0.3.2-alpha

# install basics
RUN apt-get update && apt-get install -y bash \
    curl \
    wget \
    git \
    unzip \
    gnupg \
    software-properties-common \
    ca-certificates \
    && apt-get clean

# install python 3.12
RUN add-apt-repository ppa:deadsnakes/ppa -y && \
    apt-get update && \
    apt-get install -y python3.12 python3.12-venv python3.12-dev && \
    ln -s /usr/bin/python3.12 /usr/local/bin/python3

# install Java 17.0.13 (OpenJDK)
RUN apt-get update && apt-get install -y openjdk-17-jdk

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow

# install egapx
RUN cd /opt/ && git clone https://github.com/ncbi/egapx.git
RUN python3 -m venv /opt/egapx/virtual_env && \
    /opt/egapx/virtual_env/bin/pip install --upgrade pip && \
    /opt/egapx/virtual_env/bin/pip install -r /opt/egapx/requirements.txt
# define workdir
RUN mkdir -p /opt/egapx/work
WORKDIR /opt/egapx/work
# generated config files
RUN /opt/egapx/virtual_env/bin/python3 /opt/egapx/ui/egapx.py /opt/egapx/examples/input_D_farinae_small.yaml -o /opt/egapx/example_out || true
# configure config files
RUN sed -i "s/env.GP_HOME=.*/env.GP_HOME=\"\/img\/gp\/bin\/\"/" /opt/egapx/work/egapx_config/local.config
RUN sed -i "s|env.PATH = .*|env.PATH = \"/img/gp:/img/gp/third-party/STAR/bin/Linux_x86_64/:/img/gp/third-party/diamond/:/img/gp/third-party/seqkit/bin/:/img/gp/third-party/samtools/bin/:/img/gp/third-party/miniprot/:$PATH\"|" /opt/egapx/work/egapx_config/local.config

# Add binaries EGAPx to PATH
ENV GP_HOME="/img/gp/bin/"
ENV PATH="/img/gp:/img/gp/third-party/STAR/bin/Linux_x86_64/:/img/gp/third-party/diamond/:/img/gp/third-party/seqkit/bin/:/img/gp/third-party/samtools/bin/:/img/gp/third-party/miniprot/:$PATH"

CMD ["/bin/bash"]

