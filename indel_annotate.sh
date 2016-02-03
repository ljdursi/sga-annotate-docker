#!/bin/bash 
readonly DEFREF=/reference/genome.fa.gz
readonly DEFNTHREADS=1
readonly SGABIN=/usr/local/bin/sga

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] \
    || [ ! -f "$1" ] || [ ! -f "$2" ] || [ ! -f "$3" ]
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

# clean input vcf by breaking multiallelic calls, 
# and getting rid of a sanger header line that kills sga (ends in "=")
readonly BASE=$( basename ${INPUT_VCF} )
readonly VCFBASE=${BASE%.*}
readonly CLEAN_VCF=/tmp/clean_${VCFBASE}.vcf
zcat $INPUT_VCF \
    | /deps/vcflib/bin/vcfbreakmulti \
    | grep -v "^##.*=$" \
    > ${CLEAN_VCF}

# index BAM files if necessary
if [ ! -f ${NORMAL_BAM}.bai ] 
then
    /usr/local/bin/samtools index ${NORMAL_BAM}
fi

if [ ! -f ${TUMOUR_BAM}.bai ] 
then
    /usr/local/bin/samtools index ${TUMOUR_BAM}
fi

# intermediate file, containing the output calls
# but not containing all needed new header lines
readonly BEFORE_REHEADERING_VCF=/tmp/before_headers_${VCFBASE}.vcf

${SGABIN} somatic-variant-filters \
    --annotate-only \
    --threads=$NTHREADS \
    --tumor-bam=$TUMOUR_BAM \
    --normal-bam=$NORMAL_BAM \
    --reference=$REFERENCE \
    $CLEAN_VCF > $BEFORE_REHEADERING_VCF

# output up to the start of the INFO lines
sed -n -e '1,/^##INFO/p' ${BEFORE_REHEADERING_VCF} | head -n -1
# output new header lines
cat /usr/local/share/indel.header
# output calls
sed -n -e '/^##INFO/,$p' ${BEFORE_REHEADERING_VCF}

# clean up temporary files

rm -f $CLEAN_VCF
rm -f $BEFORE_REHEADERING_VCF
