#!/bin/bash

set -e

fastqs=$1

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

mapRes_dir="${OUTPUT_DIR}/mapping_result"
mkdir -p $mapRes_dir
curr_dir=`dirname $0`

if [ $MAPPING_METHOD = "bwa" ];then
     ${curr_dir}/mapping_bwa.sh $fastqs $2 $3
elif [ $MAPPING_METHOD = "bowtie" ];then
     ${curr_dir}/mapping_bowtie.sh $fastqs $2 $3
else
     ${curr_dir}/mapping_bowtie2.sh $fastqs $2 $3
fi

## sort
echo "Sorting bam file"

ncore=$(nproc --all)
ncore=$(($ncore - 1))
mkdir -p ${mapRes_dir}/tmp
${SAMTOOLS_PATH}/samtools sort -T ${mapRes_dir}/tmp/ -@ $ncore -n -o ${mapRes_dir}/${OUTPUT_PREFIX}.sorted.bam ${mapRes_dir}/${OUTPUT_PREFIX}.bam
rm ${mapRes_dir}/${OUTPUT_PREFIX}.bam

## to mark duplicates
${SAMTOOLS_PATH}/samtools fixmate -@ $ncore -m ${mapRes_dir}/${OUTPUT_PREFIX}.sorted.bam ${mapRes_dir}/${OUTPUT_PREFIX}.fixmate.bam
rm ${mapRes_dir}/${OUTPUT_PREFIX}.sorted.bam

# Markdup needs position order
${SAMTOOLS_PATH}/samtools sort -@ $ncore -T ${mapRes_dir}/tmp/ -o ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort0.bam ${mapRes_dir}/${OUTPUT_PREFIX}.fixmate.bam
rm ${mapRes_dir}/${OUTPUT_PREFIX}.fixmate.bam

## mark duplicates
${SAMTOOLS_PATH}/samtools markdup -@ $ncore ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort0.bam ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.bam
rm ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort0.bam


${SAMTOOLS_PATH}/samtools index -@ $ncore ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.bam 



## filtering low quality and/or deplicates for downstreame analysis
${SAMTOOLS_PATH}/samtools view -f 0x2 -b -h -q 30 -@ $ncore ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.bam -o ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.MAPQ30.bam 
${SAMTOOLS_PATH}/samtools index -@ $ncore ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.MAPQ30.bam 


if [ $MAPQ -ne 30 ]; then
     ${SAMTOOLS_PATH}/samtools view -f 0x2 -b -h -q $MAPQ -@ $ncore ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.bam -o ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.MAPQ${MAPQ}.bam 
fi


## mapping stats
echo "Summarizing mapping stats ..."

curr_dir=`dirname $0`
qc_dir=${OUTPUT_DIR}/summary
mkdir -p $qc_dir
bash ${curr_dir}/mapping_qc.sh ${mapRes_dir}  $2 $3

echo "Simple mapping stats summary Done!"

echo "Getting txt file for read pair (fragment) information"
${PERL_PATH}/perl ${curr_dir}/src/simply_bam2frags.pl --read_file ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.MAPQ${MAPQ}.bam \
        --output_file ${qc_dir}/${OUTPUT_PREFIX}.fragments.tsv --samtools_path $SAMTOOLS_PATH

${TABIX_PATH}/bgzip -f ${qc_dir}/${OUTPUT_PREFIX}.fragments.tsv.len 

## sort and index fragment file
# faster sort by R data.table
${R_PATH}/Rscript --vanilla ${curr_dir}/src/sort_frags.R ${qc_dir}/${OUTPUT_PREFIX}.fragments.tsv

## index fragment file
#sort -k1,1 -k2,2n -T ${mapRes_dir}/tmp/  ${qc_dir}/${OUTPUT_PREFIX}.fragments.tsv > ${qc_dir}/${OUTPUT_PREFIX}.fragments.sorted.tsv
#mv ${qc_dir}/${OUTPUT_PREFIX}.fragments.sorted.tsv ${qc_dir}/${OUTPUT_PREFIX}.fragments.tsv
${TABIX_PATH}/bgzip  ${qc_dir}/${OUTPUT_PREFIX}.fragments.tsv 
${TABIX_PATH}/tabix -p bed ${qc_dir}/${OUTPUT_PREFIX}.fragments.tsv.gz


