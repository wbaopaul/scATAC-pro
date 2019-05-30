#!/bin/bash

input_bam=$1

filterBam_dir="${2}/filtered_bam"
mkdir -p $filterBam_dir

ncore=$(nproc --all)
ncore=$(($ncore - 1))

## rm duplicates
${SAMTOOLS_PATH}/samtools markdup -@ $ncore -r $input_bam ${filterBam_dir}/${OUTPUT_PREFIX}.dedup.bam 

## filter by MAPQ
${SAMTOOLS_PATH}/samtools view -@ $ncore -h -b -q $MAPQ ${filterBam_dir}/${OUTPUT_PREFIX}.dedup.bam -o ${filterBam_dir}/${OUTPUT_PREFIX}.dedup.MAPQ${MAPQ}.bam

echo "Filtering Bam Done!"





