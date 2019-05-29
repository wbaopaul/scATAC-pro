#!/bin/bash

set -e

input_bam=$1
qc_dir="${2}/qc_result"
mkdir -p $qc_dir


## conver bam to sam
${SAMTOOLS_PATH}/samtools view $input_bam > ${qc_dir}/tmp.sam

## use perl script to get the matrix
curr_dir=`dirname $0`

file_prefix=${OUTPUT_PREFIX}.${MAPPING_METHOD}

${PERL_PATH}/perl ${curr_dir}/cal_frac_mito.pl --read_file ${qc_dir}/tmp.sam  --output_file ${qc_dir}/${file_prefix}.freq.mito


peaks=${2}/peaks/${OUTPUT_PREFIX}_peaks.narrowPeak
${PERL_PATH}/perl ${curr_dir}/overlap_per_barcode.pl --region_file $peaks --read_file ${qc_dir}/tmp.sam --read_length 50 --output_file ${qc_dir}/${file_prefix}.overlapWith.peaks



${PERL_PATH}/perl ${curr_dir}/overlap_per_barcode.pl --region_file $PROMOTERS --read_file ${qc_dir}/tmp.sam --read_length 50 --output_file ${qc_dir}/${file_prefix}.overlapWith.promoters



${PERL_PATH}/perl ${curr_dir}/overlap_per_barcode.pl --region_file $ENHANCERS --read_file ${qc_dir}/tmp.sam --read_length 50 --output_file ${qc_dir}/${file_prefix}.overlapWith.enhancers


## merger by R
#${R_PATH}/R --vanilla --args $qc_dir $file_prefix < summary_overlap.R



#rm ${qc_dir}/tmp.sam


