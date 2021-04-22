#!/bin/bash

set -e

bam_file=$1

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

mapRes_dir="${OUTPUT_DIR}/mapping_result"
mkdir -p $mapRes_dir
curr_dir=`dirname $0`

## check whether the bam file is position sorted or not
isort=`samtools view $bam_file -H | grep HD | cut -f3`
ncore=$(nproc --all)
ncore=$(($ncore - 1))
 ## if the input is position sorted, suppose it's duplicates marked
if [[ "$isort" != "SO:coordinate"  ]]; then
    ## sort
    echo "Sorting bam file"

    mkdir -p ${mapRes_dir}/tmp
    ${SAMTOOLS_PATH}/samtools sort -T ${mapRes_dir}/tmp/ -@ $ncore -n -o ${mapRes_dir}/${OUTPUT_PREFIX}.sorted.bam $bam_file

    ## to mark duplicates
    ${SAMTOOLS_PATH}/samtools fixmate -@ $ncore -m ${mapRes_dir}/${OUTPUT_PREFIX}.sorted.bam ${mapRes_dir}/${OUTPUT_PREFIX}.fixmate.bam
    rm ${mapRes_dir}/${OUTPUT_PREFIX}.sorted.bam

    # Markdup needs position order
    ${SAMTOOLS_PATH}/samtools sort -@ $ncore -T ${mapRes_dir}/tmp/ -o ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort0.bam ${mapRes_dir}/${OUTPUT_PREFIX}.fixmate.bam
    rm ${mapRes_dir}/${OUTPUT_PREFIX}.fixmate.bam

    ## mark duplicates
    ${SAMTOOLS_PATH}/samtools markdup -@ $ncore ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort0.bam ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.bam
    rm ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort0.bam
    position_sort_bam=${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.bam
else
    position_sort_bam=$bam_file
fi

${SAMTOOLS_PATH}/samtools index -@ $ncore $position_sort_bam 


## filtering low quality and/or deplicates for downstreame analysis
${SAMTOOLS_PATH}/samtools view -f 0x2 -b -h -q 30 -@ $ncore $position_sort_bam -o ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.MAPQ30.bam 
${SAMTOOLS_PATH}/samtools index -@ $ncore ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.MAPQ30.bam 


if [ $MAPQ -ne 30 ]; then
     ${SAMTOOLS_PATH}/samtools view -f 0x2 -b -h -q $MAPQ -@ $ncore $position_sort_bam -o ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.MAPQ${MAPQ}.bam 
fi


## mapping stats
echo "Summarizing mapping stats ..."

qc_dir=${OUTPUT_DIR}/summary
mkdir -p $qc_dir
bash ${curr_dir}/mapping_qc.sh ${mapRes_dir}  $2 $3

echo "Summarize global mapping stats done!"


