#!/bin/bash



set -e
unset PYTHONPATH

input_fastqs=$1

source $2
output_dir=${OUTPUT_DIR}/demplxed_fastq
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
#dex_prefix=${fastqs[0]}
#dex_prefix=${dex_prefix##*/}
dex_prefix=$(basename ${fastqs[0]})
${PYTHON_PATH}/python ${curr_dir}/src/dex_fastq.py ${fastqs[0]} ${output_dir}/demplxed_${dex_prefix}  ${fastqs[1]}

## for the round of barcodes, add them to the read name after @, concatenate the original name by _
if [[ $kk>2 ]];then
	for (( i==2; i<=$kk; i++ ))
	do
        ${PYTHON_PATH}/python ${curr_dir}/src/dex_fastq_ul.py ${output_dir}/demplxed_${dex_prefix} ${output_dir}/demplxed_${dex_prefix}_$i ${fastqs[$i]}
        mv ${output_dir}/demplxed_${dex_prefix}_$i ${output_dir}/demplxed_${dex_prefix}
	done
fi



echo "Demultiplexing Done!"
