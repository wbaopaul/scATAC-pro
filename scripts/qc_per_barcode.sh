#!/bin/bash

set -e

input_bam=$1
qc_dir="${2}/qc_result"
mkdir -p $qc_dir


## conver bam to sam
${SAMTOOLS_PATH}/samtools view $input_bam > ${qc_dir}/tmp.sam

## use perl script to get the matrix
curr_dir=`dirname $0`

proms=/mnt/isilon/tan_lab/yuw1/scATAC-pro/annotation/hg38_promoter.bed
enhs=/mnt/isilon/tan_lab/yuw1/scATAC-pro/annotation/hg38_enhancer.bed
peaks=/mnt/isilon/tan_lab/yuw1/scATAC-pro/output/peaks/ss_peaks.narrowPeak

${PERL_PATH}/perl ${curr_dir}/cal_frac_mito.pl --read_file ${qc_dir}/tmp.sam  --output_file ${qc_dir}/freq.mito

${PERL_PATH}/perl ${curr_dir}/overlap_region_barcode.pl --region_file $peaks --read_file ${qc_dir}/tmp.sam --read_length 100 --output_file ${qc_dir}/overlapWith.peaks

${PERL_PATH}/perl ${curr_dir}/overlap_region_barcode.pl --region_file $proms --read_file ${qc_dir}/tmp.sam --read_length 100 --output_file ${qc_dir}/overlapWith.promoters

${PERL_PATH}/perl ${curr_dir}/overlap_region_barcode.pl --region_file $enhs --read_file ${qc_dir}/tmp.sam --read_length 100 --output_file ${qc_dir}/overlapWith.enhs


#rm ${qc_dir}/tmp.sam


