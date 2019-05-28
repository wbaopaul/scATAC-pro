#!/bin/bash

### will first bin genome by bedtools
set -e

input_bam=$1
mat_dir="${2}/raw_matrix"
mkdir -p $mat_dir


region_file=${mat_dir}/${SAMPLE_PREFIX}_bin.bed

${BEDTOOLS_PATH}/bedtools makewindows -g $CHROM_SIZE_FILE -w $BIN_RESL > $region_file


## conver bam to sam
${SAMTOOLS_PATH}/samtools view $input_bam > ${mat_dir}/tmp.sam

## use perl script to get the matrix
curr_dir=`dirname $0`

${PERL_PATH}/perl ${curr_dir}/get_peak_barcode_mat.pl --region_file $region_file --read_file ${mat_dir}/tmp.sam --read_length 100 --output_mat_file ${mat_dir}/${SAMPLE_PREFIX}.bin.barcode_resl${BIN_RESL}.mat  

rm ${mat_dir}/tmp.sam

echo "Get Matrix Done!"






