#!/bin/bash

input_bam=$1

filterBam_dir="${2}/filtered_bam"
mkdir -p $filterBam_dir



## rm duplicates
${SAMTOOLS_PATH}/samtools markdup -r $input_bam ${filterBam_dir}/${OUTPUT_PREFIX}.dedup.bam 

## filter by MAPQ
${SAMTOOLS_PATH}/samtools view -h -b -q $MAPQ ${filterBam_dir}/${OUTPUT_PREFIX}.dedup.bam -o ${filterBam_dir}/${OUTPUT_PREFIX}.dedup.MAPQ${MAPQ}.bam

echo "Filtering Bam Done!"





