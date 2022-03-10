#!/bin/bash



set -e
unset PYTHONPATH

input_fastqs=$1

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"


fastqs=(${input_fastqs//,/ }) ## suppose the first two fastqs is the read file, the others are index fastq files
nfile=${#fastqs[@]}
kk=$(( $nfile ))
curr_dir=`dirname $0`
mkdir -p ${OUTPUT_DIR}/demplxed_fastq
## whether the input is an 10x folder
if [[ $kk -eq 1  ]]; then
    ## cat R1, R2 and R3 for each lane
     file1=${OUTPUT_DIR}/demplxed_fastq/pooled.R1.fastq.gz
     file2=${OUTPUT_DIR}/demplxed_fastq/pooled.R2.fastq.gz
     file3=${OUTPUT_DIR}/demplxed_fastq/pooled.R3.fastq.gz
     find $1 -name '*R1_001.fastq.gz' | sort | xargs cat > $file1
     find $1 -name '*R2_001.fastq.gz' | sort | xargs cat > $file2
     find $1 -name '*R3_001.fastq.gz' | sort | xargs cat > $file3
     
     bash ${curr_dir}/dex_fastq_single.sh ${file1},${file3},${file2} $2 $3 
     
     ## cat demx fastqs from different lanes
     mv  ${OUTPUT_DIR}/demplxed_fastq/demplxed_pooled.R1.fastq.gz ${OUTPUT_DIR}/demplxed_fastq/${OUTPUT_PREFIX}.demplxed.PE1.fastq.gz
     mv  ${OUTPUT_DIR}/demplxed_fastq/demplxed_pooled.R3.fastq.gz ${OUTPUT_DIR}/demplxed_fastq/${OUTPUT_PREFIX}.demplxed.PE2.fastq.gz
    
     rm ${OUTPUT_DIR}/demplxed_fastq/pooled*
     
else
    bash ${curr_dir}/dex_fastq_single.sh $1 $2 $3
    dex_prefix1=$(basename ${fastqs[0]})
    dex_prefix2=$(basename ${fastqs[1]})
    mv ${OUTPUT_DIR}/demplxed_fastq/demplxed_${dex_prefix1} ${OUTPUT_DIR}/demplxed_fastq/${OUTPUT_PREFIX}.demplxed.PE1.fastq.gz
    mv ${OUTPUT_DIR}/demplxed_fastq/demplxed_${dex_prefix2} ${OUTPUT_DIR}/demplxed_fastq/${OUTPUT_PREFIX}.demplxed.PE2.fastq.gz
fi

echo "Demultiplexing Done!"

