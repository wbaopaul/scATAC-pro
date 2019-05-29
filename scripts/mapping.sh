#!/bin/bash

set -e


## align

mapRes_dir="${2}/mapping_result"

curr_dir=`dirname $0`

if [ $MAPPING_METHOD == "bwa" ];then
    bash ${curr_dir}/mapping_bwa.sh $1 $mapRes_dir
elif [ $MAPPING_METHOD == "bowtie" ];then
    bash ${curr_dir}/mapping_bowtie.sh $1 $mapRes_dir
else
    bash ${curr_dir}/mapping_bowtie2.sh $1 $mapRes_dir
fi




## sort
echo "Sorting bam file"
mkdir -p ${mapRes_dir}/tmp
${SAMTOOLS_PATH}/samtools sort -T ${mapRes_dir}/tmp/ -n -o ${mapRes_dir}/${SAMPLE_PREFIX}.sorted.bam ${mapRes_dir}/${SAMPLE_PREFIX}.${MAPPING_METHOD}.bam


## to mark duplicates
${SAMTOOLS_PATH}/samtools fixmate -m ${mapRes_dir}/${SAMPLE_PREFIX}.sorted.bam ${mapRes_dir}/${SAMPLE_PREFIX}.fixmate.bam


# Markdup needs position order
${SAMTOOLS_PATH}/samtools sort -o ${mapRes_dir}/${SAMPLE_PREFIX}.${MAPPING_METHOD}.positionsort.bam ${mapRes_dir}/${SAMPLE_PREFIX}.fixmate.bam



## mark duplicates
${SAMTOOLS_PATH}/samtools markdup ${mapRes_dir}/${SAMPLE_PREFIX}.${MAPPING_METHOD}.positionsort.bam ${mapRes_dir}/${SAMPLE_PREFIX}.${MAPPING_METHOD}.markdup.bam


${SAMTOOLS_PATH}/samtools index ${mapRes_dir}/${SAMPLE_PREFIX}.${MAPPING_METHOD}.markdup.bam 

## mapping stats
echo "Summarizing mapping stats ..."

curr_dir=`dirname $0`
qc_dir=${2}/qc_result
mkdir -p $qc_dir
echo ${mapRes_dir}/${SAMPLE_PREFIX}.${MAPPING_METHOD}.markdup.bam
bash $curr_dir/mapping_qc.sh ${mapRes_dir}/${SAMPLE_PREFIX}.${MAPPING_METHOD}.markdup.bam  $qc_dir



rm ${mapRes_dir}/${SAMPLE_PREFIX}.sorted.bam
rm ${mapRes_dir}/${SAMPLE_PREFIX}.fixmate.bam
rm ${mapRes_dir}/${SAMPLE_PREFIX}.${MAPPING_METHOD}.bam


echo "Simple mapping stats summary Done!"



