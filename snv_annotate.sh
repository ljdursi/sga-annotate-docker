#!/bin/bash 
readonly DEFREF=/reference/genome.fa.gz
readonly DEFNTHREADS=1
readonly READCOUNTSBIN=/usr/bin/bam-readcounts

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] \
    || [ ! -f "$1" ] || [ ! -f "$2" ] || [! -f "$3" ]
then
    echo "$0 - performs call annotation of SNVs based on bam files;"
    echo "     outputs to stdout"
    echo "Usage: $0 INPUT.VCF NORMAL.BAM TUMOUR.BAM [REFERENCE=$DEFREF]"
    echo "invocation was: $0 $1 $2 $3 $4"
    exit 
fi

readonly INPUT_VCF=$1
readonly NORMAL_BAM=$2
readonly TUMOUR_BAM=$3
readonly REFERENCE=${4-$DEFREF}

# output readcounts from input vcf

readonly NORMAL_READCOUNTS=${INPUT_VCF/.vcf/.normal.rc}
readonly TUMOUR_READCOUNTS=${INPUT_VCF/.vcf/.tumour.rc}

/usr/bin/bam-readcount --reference-fasta ${REFERENCE} --site-list <( zcat ${INPUT_VCF} | awk '{ printf "%s\t%d\t%d\n",$1,$2,$2+1 }' ) \
    --max-count 8000 --max-warnings 0 $NORMAL_BAM > ${NORMAL_READCOUNTS}

/usr/bin/bam-readcount --reference-fasta ${REFERENCE} --site-list <( zcat ${INPUT_VCF} | awk '{ printf "%s\t%d\t%d\n",$1,$2,$2+1 }' ) \
    --max-count 8000 --max-warnings 0 $TUMOUR_BAM > ${TUMOUR_READCOUNTS}

/usr/local/bin/annotate_from_readcounts.py ${INPUT_VCF} ${NORMAL_READCOUNTS} ${TUMOUR_READCOUNTS}
