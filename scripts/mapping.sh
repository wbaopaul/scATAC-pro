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

if [ $MAPPING_METHOD == "bwa" ];then
     ${curr_dir}/mapping_bwa.sh $fastqs $mapRes_dir
elif [ $MAPPING_METHOD == "bowtie" ];then
     ${curr_dir}/mapping_bowtie.sh $fastqs $mapRes_dir
else
     ${curr_dir}/mapping_bowtie2.sh $fastqs $mapRes_dir
fi


#echo "Removing reads in blacklist ..."    ## do it use fragment file, will be much faster
#${BEDTOOLS_PATH}/intersect -v -a ${mapRes_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}.bam -b $BLACKLIST > ${mapRes_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}.rmBlacklist.bam
#mv ${mapRes_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}.rmBlacklist.bam ${mapRes_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}.bam

## sort
echo "Sorting bam file"

ncore=$(nproc --all)
ncore=$(($ncore - 1))
mkdir -p ${mapRes_dir}/tmp
${SAMTOOLS_PATH}/samtools sort -T ${mapRes_dir}/tmp/ -@ $ncore -n -o ${mapRes_dir}/${OUTPUT_PREFIX}.sorted.bam ${mapRes_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}.bam


## to mark duplicates
${SAMTOOLS_PATH}/samtools fixmate -@ $ncore -m ${mapRes_dir}/${OUTPUT_PREFIX}.sorted.bam ${mapRes_dir}/${OUTPUT_PREFIX}.fixmate.bam
rm ${mapRes_dir}/${OUTPUT_PREFIX}.sorted.bam

# Markdup needs position order
${SAMTOOLS_PATH}/samtools sort -@ $ncore -T ${mapRes_dir}/tmp/ -o ${mapRes_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}.positionsort0.bam ${mapRes_dir}/${OUTPUT_PREFIX}.fixmate.bam
rm ${mapRes_dir}/${OUTPUT_PREFIX}.fixmate.bam

## mark duplicates
${SAMTOOLS_PATH}/samtools markdup -@ $ncore ${mapRes_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}.positionsort0.bam ${mapRes_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}.positionsort.bam
rm ${mapRes_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}.positionsort0.bam


${SAMTOOLS_PATH}/samtools index -@ $ncore ${mapRes_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}.positionsort.bam 



## filtering low quality and/or deplicates for downstreame analysis
${SAMTOOLS_PATH}/samtools view -f 0x2 -b -h -q 30 -@ $ncore ${mapRes_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}.positionsort.bam -o ${mapRes_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}.positionsort.MAPQ30.bam 
${SAMTOOLS_PATH}/samtools index -@ $ncore ${mapRes_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}.positionsort.MAPQ30.bam 



## mapping stats
echo "Summarizing mapping stats ..."

curr_dir=`dirname $0`
qc_dir=${OUTPUT_DIR}/summary
mkdir -p $qc_dir
bash ${curr_dir}/mapping_qc.sh ${mapRes_dir}  $qc_dir


if [ $MAPQ -ne 30 ]; then
     ${SAMTOOLS_PATH}/samtools view -f 0x2 -b -h -q $MAPQ -@ $ncore ${mapRes_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}.positionsort.bam -o ${mapRes_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}.positionsort.MAPQ${MAPQ}.bam 
fi



echo "Simple mapping stats summary Done!"

echo "Getting bed file for read pair (fragment) information"
${PERL_PATH}/perl ${curr_dir}/src/simply_bam2frags.pl --read_file ${mapRes_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}.positionsort.MAPQ${MAPQ}.bam \
        --output_file ${qc_dir}/fragments.bed --samtools_path $SAMTOOLS_PATH

echo "Remove duplicates"
${SAMTOOLS_PATH}/samtools markdup -@ $ncore -r ${mapRes_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}.positionsort.MAPQ30.bam ${mapRes_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}.positionsort.MAPQ30.noDuplicates.bam 
