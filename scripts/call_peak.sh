#!/bin/bash

input_bam=$1

peaks_dir="${2}/peaks"
mkdir -p $peaks_dir



## rm duplicates
#${MACS2_PATH}/macs2 callpeak -t $input_bam --outdir $peaks_dir -f BAM $MACS2_OPTS
${MACS2_PATH}/macs2 callpeak -t $input_bam --outdir $peaks_dir -f BAM $MACS2_OPTS --nomodel --extsize 147
## filter by MAPQ





