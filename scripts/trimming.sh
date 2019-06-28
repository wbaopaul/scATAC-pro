#!/bin/bash

set -e

input_fastqs=$1

source $2


output_dir=${OUTPUT_DIR}/trimmed_fastq
mkdir -p $output_dir

fastqs=(${input_fastqs//,/ }) ## suppose the first fastq is the read file, the others are index fastq files
nfile=${#fastqs[@]}
kk=$(( $nfile ))

if [[ $kk < 2 ]];then
  echo -e "Erro: Provide at least two fastq files: with the first one as read fastq, the others as index fastq files" >&2
  exit
fi

curr_dir=`dirname $0`


## the first barcode was add to the read name after @, and : was used to concatenate to the original name
prefix0=${fastqs[0]}
prefix0=${prefix0##*/}

prefix1=$(basename ${fastqs[1]})
ncore=$(nproc --all)
ncore=$(($nproc - 1))

if [[ $ADAPTER_SEQ != '' ]];then
    java -jar ${Trimmomatic_PATH}/*jar PE -threads 4 ${fastqs[0]} ${fastqs[1]} \
        ${output_dir}/trimmed_paired_${prefix0} ${output_dir}/trimmed_unpaired_${prefix0} \
    ${output_dir}/trimmed_paired_${prefix1} ${output_dir}/trimmed_unpaired_${prefix1} \
    ILLUMINACLIP:${ADAPTER_SEQ}:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:25
else
    echo "ADAPTER sequence not detected, using trim_galore ..."
    ${Trimgalore_PATH}/trim_galore -j 4 -o $output_dir  ${fastqs[0]} ${fastqs[1]}   --paired
fi
    
echo "Trimming Done!"
