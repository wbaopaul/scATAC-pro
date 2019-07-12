#!/bin/bash


## align
input_ff=$1
fastqs=(${input_ff//,/ })

mapRes_dir="${2}"

if [[ -z "$BWA_PATH" || ! -d "$BWA_PATH" ]];then
    echo "$BWA_PATH not found: please check you bwa installation path"
    exit
fi

## index

#echo "Indexing genome ... "
#${BWA_PATH}/bwa index $BWA_INDEX

echo "Starting alignment ... "
if [[ -z "$BWA_INDEX" || ! -f "$BWA_INDEX" ]];then
    echo "$BWA_INDEX not found: please check you bwa index file path"
    exit
fi
${BWA_PATH}/bwa mem $BWA_INDEX $BWA_OPTS ${fastqs[0]} ${fastqs[1]}  > ${mapRes_dir}/${OUTPUT_PREFIX}.bwa.sam

echo "BWA Mapping Done!"


## convert to bam
ncore=$(nproc --all)
ncore=$(($ncore - 1))

echo "Converting sam to bam ... "
${SAMTOOLS_PATH}/samtools view -@ $ncore -h -bS ${mapRes_dir}/${OUTPUT_PREFIX}.bwa.sam > ${mapRes_dir}/${OUTPUT_PREFIX}.bwa.bam

rm ${mapRes_dir}/${OUTPUT_PREFIX}.bwa.sam
