#!/bin/bash

## extract bam files for given barcodes
## the input should include a original bam file and a barcode file, seprated by comma


output_dir=${2}/extract_bam_gbarcodes
mkdir -p $output_dir

inputs=(${1//,/ })

curr_dir=`dirname $0`

input_bam=${inputs[0]}
input_barcodes=${inputs[1]}


${SAMTOOLS_PATH}/samtools view $input_bam > ${output_dir}/tmp.sam

${PERL_PATH}/perl ${curr_dir}/extract_reads_gbarcodes.pl --read_file ${output_dir}/tmp.sam --region_file $input_barcodes --output_file ${output_dir}/sam4${input_barcodes}.sam


${SAMTOOLS_PATH}/samtools view -H $input_bam > ${output_dir}/tmp.header

cat ${output_dir}/tmp.header ${output_dir}/sam4${input_barcodes}.sam > ${output_dir}/tmp.sam
 
${SAMTOOLS_PATH}/samtools view -h ${output_dir}/tmp.sam > ${output_dir}/${input_barcodes}.bam



rm ${output_dir}/*sam


