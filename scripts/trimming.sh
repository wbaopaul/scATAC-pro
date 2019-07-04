#!/bin/bash
set -e

curr_dir=`dirname $0`
#curr_dir=`cd "$curr_dir" ; pwd`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

input_fastqs=$1



output_dir=${OUTPUT_DIR}/trimmed_fastq
mkdir -p $output_dir

fastqs=(${input_fastqs//,/ }) ## suppose the first fastq is the read file, the others are index fastq files
nfile=${#fastqs[@]}
kk=$(( $nfile ))

if [[ $kk < 2 ]];then
  echo -e "Erro: Provide at least two fastq files: with the first one as read fastq, the others as index fastq files" >&2
  exit
fi



## the first barcode was add to the read name after @, and : was used to concatenate to the original name
prefix0=${fastqs[0]}
prefix0=${prefix0##*/}

prefix1=$(basename ${fastqs[1]})
ncore=$(nproc --all)
ncore=$(($nproc - 1))

if [ ${TRIM_METHOD} = 'Trimmomatic' ]; then
    echo "Using Trimmomatic ..."  
    if [ -z $TRIMMOMATIC_PATH ]; then 
        echo "error: no trimmomatic_path found" >&2
        exit 
    fi
    java -jar ${TRIMMOMATIC_PATH}/*jar PE -threads 4 ${fastqs[0]} ${fastqs[1]} \
        ${output_dir}/trimmed_paired_${prefix0} ${output_dir}/trimmed_unpaired_${prefix0} \
    ${output_dir}/trimmed_paired_${prefix1} ${output_dir}/trimmed_unpaired_${prefix1} \
    ILLUMINACLIP:${ADAPTER_SEQ}:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:25
else
    echo "Using trim_galore ..." 
    ${TRIM_GALORE_PATH}/trim_galore -j 4 -o $output_dir  ${fastqs[0]} ${fastqs[1]}   --paired
fi
    
echo "Trimming Done!" 


