#!/bin/bash 
readonly DEFREF=/reference/genome.fa.gz
readonly DEFNTHREADS=1
readonly READCOUNTSBIN=/usr/bin/bam-readcounts

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] \
    || [ ! -f "$1" ] || [ ! -f "$2" ] || [ ! -f "$3" ]
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

# clean input vcf by breaking multiallelic calls, 
# and getting rid of an sanger header line that kills some tools (ends in "=")
readonly BASE=$( basename ${INPUT_VCF} )
readonly VCFBASE=${BASE%.*}
readonly CLEAN_VCF=/tmp/clean_${VCFBASE}.vcf
zcat $INPUT_VCF \
    | /deps/vcflib/bin/vcfbreakmulti \
    | grep -v "^##.*=$" \
    > ${CLEAN_VCF}

# generate 1-based regions files from the cleaned inputs
readonly REGIONS=/tmp/regions_${VCFBASE}.txt
awk '$1 !~ /^#/{ printf "%s\t%d\t%d\n",$1,$2,$2+1 }' $CLEAN_VCF > $REGIONS

# index BAM files if necessary
if [ ! -f ${NORMAL_BAM}.bai ] 
then
    /usr/local/bin/samtools index ${NORMAL_BAM}
fi

if [ ! -f ${TUMOUR_BAM}.bai ] 
then
    /usr/local/bin/samtools index ${TUMOUR_BAM}
fi

# generate readcounts from input vcf
readonly NORMAL_READCOUNTS=/tmp/${VCFBASE}.normal.rc
readonly TUMOUR_READCOUNTS=/tmp/${VCFBASE}.tumour.rc

/usr/bin/bam-readcount --reference-fasta ${REFERENCE} \
    --site-list $REGIONS \
    --max-count 8000 $NORMAL_BAM > ${NORMAL_READCOUNTS} 2> /dev/null

/usr/bin/bam-readcount --reference-fasta ${REFERENCE} \
    --site-list $REGIONS \
    --max-count 8000 $TUMOUR_BAM > ${TUMOUR_READCOUNTS} 2> /dev/null

# intermediate file, containing the output calls
# but not containing all needed new header lines
readonly BEFORE_REHEADERING_VCF=/tmp/before_headers_${VCFBASE}.vcf

/usr/local/bin/annotate_from_readcounts.py \
    ${CLEAN_VCF} \
    ${NORMAL_READCOUNTS} ${TUMOUR_READCOUNTS} \
    > ${BEFORE_REHEADERING_VCF}

# output up to the start of the INFO lines
sed -n -e '1,/^##INFO/p' ${BEFORE_REHEADERING_VCF} | head -n -1
# output new header lines
cat /usr/local/share/snv.header
# output calls
sed -n -e '/^##INFO/,$p' ${BEFORE_REHEADERING_VCF}

# delete intermediate files
rm ${NORMAL_READCOUNTS}
rm ${TUMOUR_READCOUNTS}
rm ${REGIONS}
rm ${CLEAN_VCF}
rm ${BEFORE_REHEADERING_VCF}
