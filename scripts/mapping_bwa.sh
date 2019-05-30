#!/bin/bash

## index

#echo "Indexing genome ... "
#${MAPPING_PATH}/bwa index $BWA_INDEX


## align
input_ff=$1
fastqs=(${input_ff//,/ })

mapRes_dir="${2}"


echo "Starting alignment ... "
${MAPPING_PATH}/bwa mem $BWA_INDEX $BWA_OPTS ${fastqs[0]} ${fastqs[1]}  > ${mapRes_dir}/${OUTPUT_PREFIX}.bwa.sam

echo "BWA Mapping Done!"


## convert to bam
ncore=$(nproc --all)
ncore=$(($ncore - 1))

echo "Converting sam to bam ... "
${SAMTOOLS_PATH}/samtools view -@ $ncore -h -bS ${mapRes_dir}/${OUTPUT_PREFIX}.bwa.sam > ${mapRes_dir}/${OUTPUT_PREFIX}.bwa.bam

rm ${mapRes_dir}/${OUTPUT_PREFIX}.bwa.sam
