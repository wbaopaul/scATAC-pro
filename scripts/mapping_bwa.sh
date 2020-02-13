#!/bin/bash


## align
input_ff=$1
fastqs=(${input_ff//,/ })


# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"
isSingleEnd=${isSingleEnd^^}
mapRes_dir="${OUTPUT_DIR}/mapping_result"

if [[ -z "$BWA_PATH" || ! -d "$BWA_PATH" ]];then
    echo "$BWA_PATH not found: please check you bwa installation path"
    exit
fi


echo "Starting alignment ... "
if [[ -z "$BWA_INDEX" || ! -f "$BWA_INDEX" ]];then
    echo "Fasta file not found: please specifie the .fa file path in the user configure file"
    exit
fi

BWA_INDEX_PATH=$(dirname ${BWA_INDEX})
BWA_FA=$(basename ${BWA_INDEX})
BWA_AMB=$BWA_INDEX_PATH/${BWA_FA}.amb

if [[ ! -f "$BWA_AMB" ]];then
    echo "No index file detected: I will conduct bwa index, it will takes a few hours..."
    echo "Indexing genome ... "
    ${BWA_PATH}/bwa index $BWA_INDEX
fi


if [[ "$isSingleEnd" = "TRUE" ]]; then
    ${BWA_PATH}/bwa mem $BWA_INDEX $BWA_OPTS ${fastqs[0]}  > ${mapRes_dir}/${OUTPUT_PREFIX}.sam
else
    ${BWA_PATH}/bwa mem $BWA_INDEX $BWA_OPTS ${fastqs[0]} ${fastqs[1]}  > ${mapRes_dir}/${OUTPUT_PREFIX}.sam
fi

echo "BWA Mapping Done!"


echo "convert to bam ... "
ncore=$(nproc --all)
ncore=$((${ncore}/2))

echo "Converting sam to bam ... "
${SAMTOOLS_PATH}/samtools view -@ $ncore -h -bS ${mapRes_dir}/${OUTPUT_PREFIX}.sam > ${mapRes_dir}/${OUTPUT_PREFIX}.bam

rm ${mapRes_dir}/${OUTPUT_PREFIX}.sam
