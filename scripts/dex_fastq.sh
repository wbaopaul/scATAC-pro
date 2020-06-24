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

## whether the input is an 10x folder
if [[ $kk -eq 1  ]]; then
    ## construct R1, R2 and R3 for each lane
     R1_fastqs=($(find $1 -name '*R1*.fastq.gz' | sort ))
     for file1 in "${R1_fastqs[@]}" 
     do
        file2=${file1/R1/R2}
        file3=${file1/R1/R3}
        bash ${curr_dir}/dex_fastq_single.sh ${file1},${file3},${file2} $2 $3 &
     done
     wait
     ## cat demx fastqs from different lanes
     find ${OUTPUT_DIR}/demplxed_fastq -name '*R1*.fastq.gz' | sort | xargs cat > ${OUTPUT_DIR}/demplxed_fastq/${OUTPUT_PREFIX}.demplxed.PE1.fastq.gz
     find ${OUTPUT_DIR}/demplxed_fastq -name '*R3*.fastq.gz' | sort | xargs cat > ${OUTPUT_DIR}/demplxed_fastq/${OUTPUT_PREFIX}.demplxed.PE2.fastq.gz
     rm ${OUTPUT_DIR}/demplxed_fastq/demplxed*
else
    bash ${curr_dir}/dex_fastq_single.sh $1 $2 $3
    dex_prefix1=$(basename ${fastqs[0]})
    dex_prefix2=$(basename ${fastqs[1]})
    mv ${OUTPUT_DIR}/demplxed_fastq/demplxed_${dex_prefix1} ${OUTPUT_DIR}/demplxed_fastq/${OUTPUT_PREFIX}.demplxed.PE1.fastq.gz
    mv ${OUTPUT_DIR}/demplxed_fastq/demplxed_${dex_prefix2} ${OUTPUT_DIR}/demplxed_fastq/${OUTPUT_PREFIX}.demplxed.PE2.fastq.gz
fi

echo "Demultiplexing Done!"

