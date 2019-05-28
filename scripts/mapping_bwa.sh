#!/bin/bash

## index

#echo "Indexing genome ... "
#${MAPPING_PATH}/bwa index $GENOME_FASTA


## align
input_ff=$1
fastqs=(${input_ff//,/ })

mapRes_dir="${2}"


echo "Starting alignment ... "
${MAPPING_PATH}/bwa mem $GENOME_FASTA ${fastqs[0]} ${fastqs[1]} $MAPPING_OPTS > ${mapRes_dir}/${SAMPLE_PREFIX}.bwa.sam

echo "BWA Mapping Done!"


## convert to bam
echo "Converting sam to bam ... "
${SAMTOOLS_PATH}/samtools view -h -bS ${mapRes_dir}/${SAMPLE_PREFIX}.bwa.sam > ${mapRes_dir}/${SAMPLE_PREFIX}.bwa.bam

rm ${mapRes_dir}/${SAMPLE_PREFIX}.bwa.sam
