
############################################################
# Dockerfile for scATAC-pro pipeline
# Based on Ubuntu 20.04
############################################################

# Set the base image to Ubuntu 18.04
FROM ubuntu:20.04

MAINTAINER Wenbao Yu
ARG WHEN

# To prevent time zone prompt
ENV DEBIAN_FRONTEND=noninteractive

# Install softwares from apt repo
RUN apt-get update && apt-get install -y \
    libncurses5-dev libcurl4-openssl-dev libssl-dev libfreetype6-dev zlib1g-dev \
    git wget unzip pkg-config libboost-dev libhdf5-dev \
    default-jre apt-transport-https tabix vim pandoc \
    r-base libimage-magick-perl libgsl-dev libxml2 libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Make directory for all softwares
RUN mkdir /software
ENV PATH="/software:${PATH}"


COPY install_Rpackages_docker.R /software/install_Rpackages_docker0.R
WORKDIR /software
COPY .Renviron /software/.Renviron
RUN Rscript -e "source('install_Rpackages_docker0.R')"
RUN rm -f install_Rpackages_docker0.R

RUN apt-get update && apt-get install -y python2.7  python-pip \
    && apt-get autoremove \
    && apt-get clean
ENV PATH="/usr/bin:${PATH}"

RUN pip install --upgrade pip 
RUN pip install --upgrade cython scipy numpy 
RUN pip install --upgrade RGT 

RUN echo `which R`

## Install conda and rgt-hint
RUN wget -O tmp.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash tmp.sh -b -p /software/conda3 && rm tmp.sh

ENV PATH="/software/conda3/bin:${PATH}"

RUN conda config --add channels bioconda
    
RUN conda clean --packages -y

# Install system/math python packages (python3)
RUN conda install numpy scipy py2bit pyBigWig pysam matplotlib deeptools==3.3.1
RUN conda install cutadapt==2.5
RUN conda install macs2

RUN conda clean --packages -y

#RUN ln -s /software/conda3/bin/python /usr/bin/python

ENV PATH="/usr/bin:${PATH}"
RUN echo `which R`
##
RUN echo `R --version`

RUN conda install -y bowtie2
RUN conda install -y samtools
RUN conda install -c bioconda bedtools
RUN conda clean --packages -y

## copy standalone pandoc

COPY pandoc /usr/bin
COPY pandoc-citeproc /usr/bin

COPY configure_user.txt /software/configure_user.txt

## install scATAC-pro 
RUN git clone --branch master --single-branch https://github.com/wbaopaul/scATAC-pro.git && \
    cd scATAC-pro && \
    make configure && make install

#RUN git clone https://github.com/wbaopaul/scATAC-pro.git && cd scATAC-pro && make configure && make install && rm -rf scATAC-pro


RUN apt-get autoremove && apt-get clean


ENV PATH="/usr/local/bin/scATAC-pro_1.0.0:${PATH}"
