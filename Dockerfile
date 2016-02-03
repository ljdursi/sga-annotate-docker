# SGA Annotate variants
#
# Version 0.1
# Runs SGA ( http://github.com/jts/sga ) somatic-variant-filters --annotate-only
#  to annotate candidate indel calls with information found from the 
#  normal and tumour BAMs
#
# The invocation of SGA is:
#    sga somatic-variant-filters --annotate-only --threads=N \
#        --tumor-bam=$TUMOUR_BAM --normal-bam=$NORMAL_BAM 
#        --reference=$REFERENCE $INPUT_VCF
#
# so that we need as input the two BAM files and candidate VCF; 
# the standard pancancer reference is downloaded as part of the docker
# build process.
#
# However, some preprocessing of $INPUT_VCF must also be done; the union
# VCF from upstream may or may not have been normalized.  So we (a) run
# bcftools norm on the input VCF, and (b) merge duplicate entries 
#
FROM ubuntu:14.04
MAINTAINER Jonathan Dursi <Jonathan.Dursi@oicr.on.ca> 
LABEL Description="Runs bam-based annotation on candidate SNV, indel calls" Vendor="OICR" Version="0.1"

VOLUME /data
WORKDIR /data

# get ubuntu packages
RUN apt-get update && \
    apt-get install -y \
        automake \
        autotools-dev \
        build-essential \
        cmake \
        git \
        libhts-dev \
        libhts0 \
        libjemalloc-dev \
        libncurses5-dev \
        libsparsehash-dev \
        libz-dev \
        python \
        python-dev \
        python-pip \
        tabix \
        wget \
        zlib1g-dev 

# build remaining dependencies:
# bamtools - for SGA
RUN mkdir -p /deps && \
    cd /deps && \
    wget https://github.com/pezmaster31/bamtools/archive/v2.4.0.tar.gz && \
    tar -xzvf v2.4.0.tar.gz && \
    rm -rf v2.4.0.tar.gz && \
    cd bamtools-2.4.0 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make

# bam-readcount: for SNV annotation
RUN mkdir -p /deps && \
    cd /deps && \
    wget https://github.com/genome/bam-readcount/archive/v0.7.4.tar.gz && \
    tar -xzvf v0.7.4.tar.gz && \
    rm v0.7.4.tar.gz && \
    cd bam-readcount-0.7.4 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install

# samtools - for indexing reference, etc
RUN mkdir -p /deps && \
    cd /deps && \
    wget https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2 && \
    tar -xjvf samtools-1.3.tar.bz2 && \
    rm samtools-1.3.tar.bz2 && \
    cd samtools-1.3 && \
    make prefix=/usr/local/ install && \
    cd .. && \
    rm -rf samtools-1.3

# vcflib - for tools like vcfbreakmulti
# set a fixed version for reproducibility
RUN mkdir -p /deps && \
    cd /deps && \
    git clone --recursive git://github.com/ekg/vcflib.git && \
    cd vcflib && \
    git checkout d453d91592fe8a74d92b49cd6c7cd73f79a8b70b && \
    make 

# get pancan standard reference
RUN mkdir -p /reference && \
    cd /reference && \
    wget ftp://ftp.sanger.ac.uk/pub/project/PanCancer/genome.fa.gz && \
    gunzip genome.fa.gz ; \
    bgzip genome.fa && \
    samtools faidx genome.fa.gz

# build SGA
RUN mkdir -p /src && \
    cd /src && \
    wget https://github.com/jts/sga/archive/v0.10.14.tar.gz && \
    tar -xzvf v0.10.14.tar.gz && \
    cd sga-0.10.14/src && \
    ./autogen.sh && \
    ./configure --with-bamtools=/deps/bamtools-2.4.0 --with-jemalloc=/usr --prefix=/usr/local && \
    make && \
    make install

# get pyvcf for annotate_from_readcounts.py
RUN pip install --upgrade pip && \
    pip install pyvcf

# Put auxilliary data in /usr/local/share
COPY indel.header /usr/local/share
COPY snv.header /usr/local/share

# Put scripts in /usr/local/bin
COPY indel_annotate.sh /usr/local/bin
COPY snv_annotate.sh /usr/local/bin
COPY annotate.sh /usr/local/bin
COPY annotate_from_readcounts.py /usr/local/bin

ENTRYPOINT ["/usr/local/bin/annotate.sh"]
