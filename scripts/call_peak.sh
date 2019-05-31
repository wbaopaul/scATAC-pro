#!/bin/bash

input_bam=$1

peaks_dir="${2}/peaks"
mkdir -p $peaks_dir

out_prefix=${OUTPUT_PREFIX}.${MAPPING_METHOD}

## call peaks
${MACS2_PATH}/macs2 callpeak -t $input_bam --outdir $peaks_dir -n $out_prefix -f BAM $MACS2_OPTS 
#${MACS2_PATH}/macs2 callpeak -t $input_bam --outdir $peaks_dir -f BAM $MACS2_OPTS --nomodel --extsize 147


## remove peaks overlapped with blacklist
${BEDTOOLS_PATH}/bedtools intersect -a ${peaks_dir}/${out_prefix}_peaks.narrowPeak -b $BLACKLIST -v \
    > ${peaks_dir}/${out_prefix}_peaks_filterBlacklist.narrowPeak





