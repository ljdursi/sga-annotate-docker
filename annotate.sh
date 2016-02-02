#!/bin/bash 
if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] \
    || [ ! -f "$2" ] || [ ! -f "$3" ] || [ ! -f "$4" ]
then
    echo "$0 - performs call annotation based on bam files;"
    echo "     outputs to stdout"
    echo "Usage: $0 {SNV|indel} INPUT.VCF NORMAL.BAM TUMOUR.BAM [options]"
    echo "invocation was: $0 $1 $2 $3 $4"
    exit 
fi


if [ ${1:0:3} == "SNV" ] || [ ${1:0:3} == "snv" ]
then
    >&2 echo "#SNV Annotation"
    /usr/local/bin/snv_annotate.sh $2 $3 $4 $5 $6
else
    >&2 echo "#indel Annotation"
    /usr/local/bin/indel_annotate.sh $2 $3 $4 $5 $6
fi
