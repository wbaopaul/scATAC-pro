#!/bin/bash

## output mapping qc results

input_bam=$1
output_dir=$2
ncore=$(nproc --all)
ncore=$(($ncore - 1))
${SAMTOOLS_PATH}/samtools flagstat $input_bam > ${2}/qc_${MAPPING_METHOD}_flagstat.txt
${SAMTOOLS_PATH}/samtools idxstats $input_bam > ${2}/qc_${MAPPING_METHOD}_idxstat.txt

tmp_sam=${2}/tmp.sam
${SAMTOOLS_PATH}/samtools view -@ $ncore -q 1 $input_bam > ${2}/tmp.sam

if [ $MAPPING_METHOD == 'bwa' ]; then
  total_uniq_map=$(grep -v -e 'XA:Z' -e 'SA:Z:' $tmp_sam | wc -l)  ## number of unique mapped reads
elif [ $MAPPING_METHOD == 'bowtie' ]; then
  total_uniq_map=$( grep -E "@|NM:" $tmp_sam | grep -v "XS:" | wc -l )
else 
  total_uniq_map=$( grep -E "@|NM:" $tmp_sam | grep -v "XS:" | wc -l )
fi


total_q10=$(${SAMTOOLS_PATH}/samtools view -@ $ncore -q 10 $input_bam | wc -l)
total_q30=$(${SAMTOOLS_PATH}/samtools view -@ $ncore -q 30 $input_bam | wc -l)


echo "total_uniq_map    $total_uniq_map" > ${2}/qc_${MAPPING_METHOD}_tmp.txt 
echo "total_q10 $total_q10" >> ${2}/qc_${MAPPING_METHOD}_tmp.txt 
echo "total_q30 $total_q30" >> ${2}/qc_${MAPPING_METHOD}_tmp.txt 


rm $tmp_sam

## using R to make the result as a table


