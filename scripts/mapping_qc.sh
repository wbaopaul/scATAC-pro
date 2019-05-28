#!/bin/bash

## output mapping qc results

input_bam=$1
output_dir=$2


${SAMTOOLS_PATH}/samtools flagstat $input_bam > ${2}/qc_flagstat.txt
${SAMTOOLS_PATH}/samtools idxstats $input_bam > ${2}/qc_idxstat.txt


total_uniq_map=`${SAMTOOLS_PATH}/samtools view -q 1 $input_bam | grep -v -e 'XA:Z' -e 'SA:Z:' | wc -l`  ## number of unique mapped reads
total_q10=`${SAMTOOLS_PATH}/samtools view -q 10 $input_bam | wc -l` 
total_q30=`${SAMTOOLS_PATH}/samtools view -q 30 $input.bam | wc -l`


echo "total_uniq_map \t $total_uniq_map" > ${2}/qc_tmp.txt 
echo "total_q10 \t $total_q10" >> ${2}/qc_tmp.txt 
echo "total_q30 \t $total_q30" >> ${2}/qc_tmp.txt 

## using R to make the result as a table


