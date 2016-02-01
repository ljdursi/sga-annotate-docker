#!/bin/bash 
readonly DEFREF=/reference/genome.fa.gz
readonly DEFNTHREADS=1
readonly SGABIN=/usr/local/bin/sga

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] \
    || [ ! -f "$1" ] || [ ! -f "$2" ] || [! -f "$3" ]
then
    echo "$0 - performs call annotation of indels based on bam files;"
    echo "     outputs to stdout"
    echo "Usage: $0 INPUT.VCF NORMAL.BAM TUMOUR.BAM [NTHREADS=$DEFNTHREADS] [REFERENCE=$DEFREF]"
    echo "invocation was: $0 $1 $2 $3 $4"
    exit 
fi

readonly INPUT_VCF=$1
readonly NORMAL_BAM=$2
readonly TUMOUR_BAM=$3
readonly NTHREADS=${4-$DEFNTHREADS}
readonly REFERENCE=${5-$DEFREF}

# todo - use vcfbreakmulti to 
${SGABIN} somatic-variant-filters \
    --annotate-only \
    --threads=$NTHREADS \
    --tumor-bam=$TUMOUR_BAM \
    --normal-bam=$NORMAL_BAM \
    --variant-type=INDEL,COMPLEX \
    --reference=$REFERENCE \
    ${INPUT_VCF}
