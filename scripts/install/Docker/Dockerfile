
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
    python3 python3-setuptools python3-pip \
    git wget unzip pkg-config libboost-dev libhdf5-dev \
    default-jre apt-transport-https tabix vim \
    r-base libimage-magick-perl libgsl-dev libxml2 libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Make directory for all softwares
RUN mkdir /software
ENV PATH="/software:${PATH}"


COPY install_Rpackages.R /software/install_Rpackages.R
WORKDIR /software
RUN Rscript -e "source('install_Rpackages_docker.R')"


# Install system/math python packages (python3)
RUN pip3 install --no-cache-dir numpy scipy py2bit pyBigWig pysam matplotlib deeptools==3.3.1  
RUN pip3 install --no-cache-dir cutadapt==2.5
RUN pip3 install --no-cache-dir macs2

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN echo `which python`
RUN echo `python --version`
RUN echo `R --version`
RUN echo `which pip`

# Install bedtools 2.29.0
#RUN git clone --branch v2.29.0 --single-branch https://github.com/arq5x/bedtools2.git && \
#    cd bedtools2 && make && make install && cd ../ && rm -rf bedtools2*

# Install samtools 1.9
RUN git clone --branch 1.9 --single-branch https://github.com/samtools/samtools.git && \
    git clone --branch 1.9 --single-branch https://github.com/samtools/htslib.git && \
    cd samtools && make && make install && cd ../ && rm -rf samtools* htslib*


RUN git clone --branch master --single-branch https://github.com/wbaopaul/scATAC-pro.git && \
    cd scATAC-pro && make configure && make install

RUN rm -rf scATAC-pro

COPY configure_user.txt /software/configure_user.txt
RUN rm -f install_Rpackages_docker.R

