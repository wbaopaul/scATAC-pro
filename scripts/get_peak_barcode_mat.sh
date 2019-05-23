#!/bin/bash

### will search peaks/*narrowPeak under $2

input_bam=$1
region_file=${2}/peaks/${SAMPLE_PREFIX}_peaks.narrowPeak
mat_dir="${2}/raw_matrix"
mkdir -p $mat_dir

## get  aggregated peak counts file (optional)
#cat $region_file | cut -f1-3 > ${mat_dir}/tmp_file
#${BEDTOOLS_PATH}/bedtools coverage -counts -a ${mat_dir}/tmp_file -b ${input_bam} > ${mat_dir}/${SAMPLE_PREFIX}.peak.count
#rm ${mat_dir}/tmp_file


## conver bam to sam
${SAMTOOLS_PATH}/samtools view $input_bam > ${mat_dir}/tmp.sam

## use perl script to get the matrix
curr_dir=`dirname $0`

${PERL_PATH}/perl ${curr_dir}/get_peak_barcode_mat.pl --region_file $region_file --read_file ${mat_dir}/tmp.sam --read_length 100 --output_file ${mat_dir}/${SAMPLE_PREFIX}.peak.barcode.mat

rm ${mat_dir}/tmp.sam

echo "Get Matrix Done!"






