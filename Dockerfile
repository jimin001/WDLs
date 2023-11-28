FROM ubuntu:18.04
MAINTAINER Jimin Park, jpark621@ucsc.edu

# update and install dependencies
RUN apt-get update && \
    apt-get -y install time git make wget autoconf gcc g++ vim sudo build-essential bzip2 zlib1g-dev libbz2-dev \
        libcurl4-gnutls-dev liblzma-dev libncurses5-dev libncursesw5-dev libssl-dev

RUN mkdir -p /home/apps

### samtools
ARG SAMTOOLS_VERSION=1.17
RUN cd /home/apps && \
	wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
	tar -vxjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
	rm -rf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
	cd samtools-${SAMTOOLS_VERSION} && \
	make
ENV PATH="/home/apps/samtools-${SAMTOOLS_VERSION}:$PATH"

### htslib
ARG HTSLIB_VERSION=1.17

RUN cd /home/apps && \
	wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    tar -vxjf htslib-${HTSLIB_VERSION}.tar.bz2 && \
    rm -rf htslib-${HTSLIB_VERSION}.tar.bz2 && \
    cd htslib-${HTSLIB_VERSION} && \
    make
ENV PATH="/home/apps/htslib-${HTSLIB_VERSION}:$PATH"


### bcftools
ARG BCFTOOLS_VERSION=1.17
RUN cd /home/apps && \
	wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar -vxjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    rm -rf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION}/ && \
    make
ENV PATH="/home/apps/bcftools-${BCFTOOLS_VERSION}:$PATH"

### minimap2
ARG MINIMAP2_VERSION=2.26
RUN cd /home/apps && \
	wget https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}.tar.bz2 && \
	tar -vxjf minimap2-${MINIMAP2_VERSION}.tar.bz2 && \
	rm minimap2-${MINIMAP2_VERSION}.tar.bz2 && \
	cd minimap2-${MINIMAP2_VERSION} && \
	make
ENV PATH="/home/apps/minimap2-${MINIMAP2_VERSION}:$PATH"

### bwa-mem2
ARG BWA_MEM_VERSION=2.2.1
RUN cd /home/apps && \
	apt-get -y update && apt-get -y install bzip2 && apt-get -y install curl && \
	curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWA_MEM_VERSION}/bwa-mem2-${BWA_MEM_VERSION}_x64-linux.tar.bz2 | tar jxf - && \
	rm -rf bwa-mem2-${BWA_MEM_VERSION}_x64-linux.tar.bz2
ENV PATH="/home/apps/bwa-mem2-2.2.1_x64-linux:$PATH"




