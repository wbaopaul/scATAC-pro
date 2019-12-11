#!/bin/bash



set -e
unset PYTHONPATH

input_fastqs=$1

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

output_dir=${OUTPUT_DIR}/demplxed_fastq
mkdir -p $output_dir

fastqs=(${input_fastqs//,/ }) ## suppose the first two fastqs is the read file, the others are index fastq files
nfile=${#fastqs[@]}
kk=$(( $nfile ))
curr_dir=`dirname $0`

isSingleEnd=${isSingleEnd^^}
if [[ "$isSingleEnd" = "TRUE" ]]; then
    
    if [[ $kk < 2 ]];then
      echo -e "Erro: Provide at least two fastq files: one for read the other for index (barcode) " >&2
      exit
    fi
    ## the first barcode was add to the read name after @, and : was used to concatenate to the original name
    dex_prefix1=$(basename ${fastqs[0]})

    ${PYTHON_PATH}/python ${curr_dir}/src/dex_fastq.py ${fastqs[0]} ${output_dir}/demplxed_${dex_prefix1}  ${fastqs[1]} 


    ## for extra round of barcodes, add them to the read name after @, concatenate the original name by _
    if [[ $kk>2 ]];then
        for (( i==2; i<=$kk; i++ ))
        do
            ${PYTHON_PATH}/python ${curr_dir}/src/dex_fastq_ul.py ${output_dir}/demplxed_${dex_prefix1} ${output_dir}/demplxed_${dex_prefix1}_$i ${fastqs[$i]} 
            mv ${output_dir}/demplxed_${dex_prefix1}_$i ${output_dir}/demplxed_${dex_prefix1}
        done
    fi
else
    if [[ $kk < 3 ]];then
      echo -e "Erro: Provide at least three fastq files: the first two for paired-end reads the other for index (barcode) " >&2
      exit
    fi
    
    ## the first barcode was add to the read name after @, and : was used to concatenate to the original name
    dex_prefix1=$(basename ${fastqs[0]})
    dex_prefix2=$(basename ${fastqs[1]})



    ${PYTHON_PATH}/python ${curr_dir}/src/dex_fastq.py ${fastqs[0]} ${output_dir}/demplxed_${dex_prefix1}  ${fastqs[2]} &

    ${PYTHON_PATH}/python ${curr_dir}/src/dex_fastq.py ${fastqs[1]} ${output_dir}/demplxed_${dex_prefix2}  ${fastqs[2]} &
    wait

    ## for the round of barcodes, add them to the read name after @, concatenate the original name by _
    if [[ $kk>3 ]];then
        for (( i=3; i<=$kk; i++ ))
        do
            ${PYTHON_PATH}/python ${curr_dir}/src/dex_fastq_ul.py ${output_dir}/demplxed_${dex_prefix1} ${output_dir}/tmp${i}_demplxed_${dex_prefix1} ${fastqs[$i]} &
            ${PYTHON_PATH}/python ${curr_dir}/src/dex_fastq_ul.py ${output_dir}/demplxed_${dex_prefix2} ${output_dir}/tmp${i}_demplxed_${dex_prefix2} ${fastqs[$i]} &
            wait
            mv ${output_dir}/tmp${i}_demplxed_${dex_prefix1} ${output_dir}/demplxed_${dex_prefix1}
            mv ${output_dir}/tmp${i}_demplxed_${dex_prefix2} ${output_dir}/demplxed_${dex_prefix2}
        done
    fi



fi


echo "Demultiplexing Done!"

