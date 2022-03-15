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
mkdir -p ${OUTPUT_DIR}/demplxed_fastq

fastqs=(${input_fastqs//,/ }) ## suppose the first fastq is the read file, the others are index fastq files
nfile=${#fastqs[@]}
kk=$(( $nfile ))
#isSingleEnd=${isSingleEnd^^}
isSingleEnd=$(echo $isSingleEnd | tr a-z A-Z)

if [[ "$isSingleEnd" = "TRUE" ]]; then
    ## make the output name consistent
    dex_fastq1=${OUTPUT_DIR}/demplxed_fastq/${OUTPUT_PREFIX}.demplxed.PE1.fastq.gz
    if [[ ! -f ${dex_fastq1} ]]; then
        ln -s ${fastqs[0]} $dex_fastq1
    fi
    prefix0=$(basename $dex_fastq1)
    if [ "$TRIM_METHOD" = 'Trimmomatic' ]; then
        echo "Using Trimmomatic ..."  
        if [ -z $TRIMMOMATIC_PATH ]; then 
            echo "error: no trimmomatic_path found" >&2
            exit 
        fi

        trimmed_fastq1=${output_dir}/trimmed_paired_${prefix0}
        if [ -f "$trimmed_fastq1" ]; then
            echo -e "Trimmed fastq file $trimmed_fastq1 exist, I will skip trimming
                    reads!"
            exit
        fi

        java -jar ${TRIMMOMATIC_PATH}/*jar SE -phred33 $dex_fastq1 ${output_dir}/trimmed_${prefix0} \
             ILLUMINACLIP:${ADAPTER_SEQ}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25

        mv $trimmed_fastq1 ${OUTPUT_DIR}/trimmed_fastq/${OUTPUT_PREFIX}.trimmed.demplxed.PE1.fastq.gz
        echo "Trimming Done!" 
    elif [ "$TRIM_METHOD" = 'trim_galore' ]; then
        echo "Using trim_galore ..." 
        unset PYTHONHOME
        unset PYTHONPATH
        trimmed_fastq1=${OUTPUT_DIR}/trimmed_fastq/${OUTPUT_PREFIX}.demplxed.PE1_val_1.fq.gz
        if [ -f "$trimmed_fastq1" ]; then
            echo -e "Trimmed fastq file $trimmed_fastq1 exist, I will skip trimming
                    reads!"
            exit
        fi
        ${TRIM_GALORE_PATH}/trim_galore -j 4 -o $output_dir  $dex_fastq1 --gzip --path_to_cutadapt ${CUTADAPT_PATH}/cutadapt
 
        mv $trimmed_fastq1 ${OUTPUT_DIR}/trimmed_fastq/${OUTPUT_PREFIX}.trimmed.demplxed.PE1.fastq.gz
        echo "Trimming Done!" 
    else
        echo "You have not specify TRIM_METHOD, so I do not trim the reads"
    fi
else
    ## make the output name consistent
    dex_fastq1=${OUTPUT_DIR}/demplxed_fastq/${OUTPUT_PREFIX}.demplxed.PE1.fastq.gz
    dex_fastq2=${OUTPUT_DIR}/demplxed_fastq/${OUTPUT_PREFIX}.demplxed.PE2.fastq.gz
    if [[ ! -f ${dex_fastq1} ]]; then
        ln -s ${fastqs[0]} $dex_fastq1
        ln -s ${fastqs[1]} $dex_fastq2
    fi
    prefix0=$(basename $dex_fastq1)
    prefix1=$(basename $dex_fastq2)
    if [ "$TRIM_METHOD" = 'Trimmomatic' ]; then
        echo "Using Trimmomatic ..."  
        if [ -z $TRIMMOMATIC_PATH ]; then 
            echo "error: no trimmomatic_path found" >&2
            exit 
        fi

        trimmed_fastq1=${output_dir}/trimmed_paired_${prefix0}
        trimmed_fastq2=${output_dir}/trimmed_paired_${prefix1}
        java -jar ${TRIMMOMATIC_PATH}/*jar PE -threads 4 $dex_fastq1 $dex_fastq2 \
            ${output_dir}/trimmed_paired_${prefix0} ${output_dir}/trimmed_unpaired_${prefix0} \
        ${output_dir}/trimmed_paired_${prefix1} ${output_dir}/trimmed_unpaired_${prefix1} \
        ILLUMINACLIP:${ADAPTER_SEQ}:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:25

        mv $trimmed_fastq1 ${OUTPUT_DIR}/trimmed_fastq/${OUTPUT_PREFIX}.trimmed.demplxed.PE1.fastq.gz
        mv $trimmed_fastq2 ${OUTPUT_DIR}/trimmed_fastq/${OUTPUT_PREFIX}.trimmed.demplxed.PE2.fastq.gz
        echo "Trimming Done!" 
    elif [ "$TRIM_METHOD" = 'trim_galore' ]; then
        echo "Using trim_galore ..." 
        unset PYTHONHOME
        unset PYTHONPATH

        ${TRIM_GALORE_PATH}/trim_galore -j 4 -o $output_dir  $dex_fastq1 $dex_fastq2 --paired --gzip --path_to_cutadapt ${CUTADAPT_PATH}/cutadapt

        trimmed_fastq1=${OUTPUT_DIR}/trimmed_fastq/${OUTPUT_PREFIX}.demplxed.PE1_val_1.fq.gz
        trimmed_fastq2=${OUTPUT_DIR}/trimmed_fastq/${OUTPUT_PREFIX}.demplxed.PE2_val_2.fq.gz

        mv $trimmed_fastq1 ${OUTPUT_DIR}/trimmed_fastq/${OUTPUT_PREFIX}.trimmed.demplxed.PE1.fastq.gz
        mv $trimmed_fastq2 ${OUTPUT_DIR}/trimmed_fastq/${OUTPUT_PREFIX}.trimmed.demplxed.PE2.fastq.gz
        echo "Trimming Done!" 
    else
        echo "You have not specify TRIM_METHOD, so I do not trim the reads"
    fi

fi
