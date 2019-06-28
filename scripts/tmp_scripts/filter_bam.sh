#!/bin/bash

#### only keep proper paired, high mapq reads
#### for down stream analysis
#### with further filter barcodes with reads less than 10

input_bam=$1

filterBam_dir="${2}/filtered_bam"
mkdir -p $filterBam_dir

ncore=$(nproc --all)
ncore=$(($ncore - 1))

out_prefix=${OUTPUT_PREFIX}.${MAPPING_METHOD}

echo "filter by MAPQ${MAPQ}..."
${SAMTOOLS_PATH}/samtools view -@ $ncore -h -b -q $MAPQ $input_bam -o ${filterBam_dir}/tmp.bam

echo "only keep proper paired pairs..."
${SAMTOOLS_PATH}/samtools view -@ $ncore -f 0x2 -h -b  ${filterBam_dir}/tmp.bam -o ${filterBam_dir}/${out_prefix}.paired.MAPQ${MAPQ}.bam

rm ${filterBam_dir}/tmp.bam

echo "remove duplicates..."
${SAMTOOLS_PATH}/samtools markdup -@ $ncore -r ${filterBam_dir}/${out_prefix}.paired.MAPQ${MAPQ}.bam ${filterBam_dir}/${out_prefix}.dedup.paired.MAPQ${MAPQ}.bam 

rm ${filterBam_dir}/${out_prefix}.paired.MAPQ${MAPQ}.bam


echo "Filtering Bam Done!"





