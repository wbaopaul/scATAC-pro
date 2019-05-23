#!/bin/bash

## index

#echo "Indexing genome ... "
#${BWA_PATH}/bwa index $GENOME_FASTA


## align
input_ff=$1
fastqs=(${input_ff//,/ })

mapRes_dir="${2}/mapping_result"


echo "Starting alignment ... "
${BWA_PATH}/bwa mem $GENOME_FASTA ${fastqs[0]} ${fastqs[1]} > ${mapRes_dir}/${SAMPLE_PREFIX}.sam

echo "BWA Mapping Done!"


## convert to bam
echo "Converting sam to bam ... "
${SAMTOOLS_PATH}/samtools view -h -bS ${mapRes_dir}/${SAMPLE_PREFIX}.sam > ${mapRes_dir}/${SAMPLE_PREFIX}.bam

rm ${mapRes_dir}/${SAMPLE_PREFIX}.sam

## sort
echo "Sorting bam file"
mkdir -p ${mapRes_dir}/tmp
${SAMTOOLS_PATH}/samtools sort -T ${mapRes_dir}/tmp/ -n -o ${mapRes_dir}/${SAMPLE_PREFIX}.sorted.bam ${mapRes_dir}/${SAMPLE_PREFIX}.bam


## to mark duplicates
${SAMTOOLS_PATH}/samtools fixmate -m ${mapRes_dir}/${SAMPLE_PREFIX}.sorted.bam ${mapRes_dir}/${SAMPLE_PREFIX}.fixmate.bam


# Markdup needs position order
${SAMTOOLS_PATH}/samtools sort -o ${mapRes_dir}/${SAMPLE_PREFIX}.positionsort.bam ${mapRes_dir}/${SAMPLE_PREFIX}.fixmate.bam



## mark duplicates
${SAMTOOLS_PATH}/samtools markdup ${mapRes_dir}/${SAMPLE_PREFIX}.positionsort.bam ${mapRes_dir}/${SAMPLE_PREFIX}.markdup.bam


${SAMTOOLS_PATH}/samtools index ${mapRes_dir}/${SAMPLE_PREFIX}.markdup.bam 

## mapping flag stats
${SAMTOOLS_PATH}/samtools flagstat ${mapRes_dir}/${SAMPLE_PREFIX}.markdup.bam > ${mapRes_dir}/${SAMPLE_PREFIX}.markdup.map.flagstat.txt


${SAMTOOLS_PATH}/samtools idxstats ${mapRes_dir}/${SAMPLE_PREFIX}.markdup.bam > ${mapRes_dir}/${SAMPLE_PREFIX}.markdup.map.idxstat.txt



rm ${mapRes_dir}/${SAMPLE_PREFIX}.sorted.bam
rm ${mapRes_dir}/${SAMPLE_PREFIX}.fixmate.bam


echo "Simple mapping stats summary Done!"
