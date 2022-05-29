# Title          :meta-index
# Description    :A pipeline for automatically indexing genomes and accurately characterizing metagenome-assembled genomes with sequence bloom trees
# Author         :Fabio Cumbo (fabio.cumbo@gmail.com)

FROM ubuntu:18.04

MAINTAINER fabio.cumbo@gmail.com

# Set the working directory
WORKDIR /home

# Upgrade installed packages
RUN apt-get update && apt-get upgrade -y && apt-get clean

# Installing basic dependancies
RUN apt-get install -y \ 
        bc \
        build-essential \
        curl \
        git \
        gzip \
        make \
        python3.7 \
        python3.7-dev \
        python3.7-venv \
        python3.7-distutils \
        unzip \
        wget

# Make Python 3.7 available with venv
RUN python3.7 -m venv /venv
ENV PATH="/venv/bin:${PATH}"

# Upgrade pip using Python 3.7
RUN curl -s https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
    python get-pip.py --force-reinstall && \
    rm get-pip.py

# Install CheckM dependencies
RUN mkdir checkm-deps
# HMMER
RUN cd /home/checkm-deps && wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz && \
    tar -xvzf hmmer-3.1b2-linux-intel-x86_64.tar.gz && \
    cd hmmer-3.1b2-linux-intel-x86_64 && ./configure &&  \
    make && make install
# prodigal
RUN cd /home/checkm-deps && wget https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux && \
    mv prodigal.linux /usr/local/bin/prodigal && chmod 755 /usr/local/bin/prodigal
# pplacer
RUN cd /home/checkm-deps && wget https://github.com/matsen/pplacer/releases/download/v1.1.alpha19/pplacer-linux-v1.1.alpha19.zip && \
    unzip pplacer-linux-v1.1.alpha19.zip && \
    cp pplacer-Linux-v1.1.alpha19/pplacer /usr/local/bin/ && \
    cp pplacer-Linux-v1.1.alpha19/guppy /usr/local/bin/ && \
    cp pplacer-Linux-v1.1.alpha19/rppr /usr/local/bin/

# Install CheckM
RUN pip install numpy && \
    pip install matplotlib && \
    pip install pysam && \
    pip install checkm-genome

# Download CheckM data 
RUN cd /home && mkdir checkm-data && cd checkm-data && \
    wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz && \
    tar -xvzf checkm_data_2015_01_16.tar.gz
# Set CheckM database folder
CMD checkm data setRoot /home/checkm-data

# Install ncbitax2lin
RUN pip install ncbitax2lin

# Install miniconda
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
# Add conda channels
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# Install howdesbt and kmtricks with conda
RUN conda install -y howdesbt kmtricks

# Download meta-index
RUN mkdir -p /home/git && cd /home/git && \
    git clone https://github.com/BlankenbergLab/meta-index && \
    cd /home/git/meta-index && \
    chmod +x meta-index && \
    cd /home/git/meta-index/modules && \
    chmod +x *
# Add meta-index to PATH
ENV PATH="${PATH}:/home/git/meta-index"

WORKDIR /home
ENTRYPOINT /bin/bash