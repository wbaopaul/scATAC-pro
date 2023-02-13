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


ncore=$(nproc --all)
ncore=$(($ncore - 1))
mkdir -p ${mapRes_dir}/tmp

## sort
echo "Sorting bam file"
## sort by read name
${SAMTOOLS_PATH}/samtools sort -m 2G -T ${mapRes_dir}/tmp/ -@ $ncore -n -o ${mapRes_dir}/${OUTPUT_PREFIX}.sorted.bam ${mapRes_dir}/${OUTPUT_PREFIX}.bam
rm ${mapRes_dir}/${OUTPUT_PREFIX}.bam

## to mark duplicates
${SAMTOOLS_PATH}/samtools fixmate -@ $ncore -m ${mapRes_dir}/${OUTPUT_PREFIX}.sorted.bam ${mapRes_dir}/${OUTPUT_PREFIX}.fixmate.bam
rm ${mapRes_dir}/${OUTPUT_PREFIX}.sorted.bam


# Markdup needs position order
${SAMTOOLS_PATH}/samtools sort -m 2G -@ $ncore -T ${mapRes_dir}/tmp/ -o ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort0.bam ${mapRes_dir}/${OUTPUT_PREFIX}.fixmate.bam
rm ${mapRes_dir}/${OUTPUT_PREFIX}.fixmate.bam

## mark duplicates
${SAMTOOLS_PATH}/samtools markdup -@ $ncore ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort0.bam ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.bam
rm ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort0.bam

${SAMTOOLS_PATH}/samtools index -@ $ncore ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.bam 

if [ -z "$SHIFT_READS_IN_BAM" ]; then
    SHIFT_READS_IN_BAM=FALSE
fi

if [ ${SHIFT_READS_IN_BAM} = 'TRUE' ]; then
    ${DEEPTOOLS_PATH}/alignmentSieve --numberOfProcessors $ncore --ATACshift --bam ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.bam -o ${mapRes_dir}/${OUTPUT_PREFIX}.shifted.bam
    rm ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.bam
    ${SAMTOOLS_PATH}/samtools sort -m 2G -@ $ncore -T ${mapRes_dir}/tmp/ -o ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.bam ${mapRes_dir}/${OUTPUT_PREFIX}.shifted.bam
    rm ${mapRes_dir}/${OUTPUT_PREFIX}.shifted.bam
    ${SAMTOOLS_PATH}/samtools index -@ $ncore ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.bam 
fi

## filtering low quality for downstreame analysis
if [ ${isSingleEnd} = 'TRUE' ]; then
    ${SAMTOOLS_PATH}/samtools view -b -h -q $MAPQ -@ $ncore ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.bam -o ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.MAPQ${MAPQ}.bam 
else
     ${SAMTOOLS_PATH}/samtools view -f 0x2 -b -h -q $MAPQ -@ $ncore ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.bam -o ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.MAPQ${MAPQ}.bam 
fi

${SAMTOOLS_PATH}/samtools index -@ $ncore ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.MAPQ${MAPQ}.bam 


## mapping stats
echo "Summarizing mapping stats ..."

curr_dir=`dirname $0`
qc_dir=${OUTPUT_DIR}/summary
mkdir -p $qc_dir
bash ${curr_dir}/mapping_qc.sh ${mapRes_dir}  $2 $3

echo "Simple mapping stats summary Done!"

echo "Getting txt file for read pair (fragment) information"

bam2frag=simply_bam2frags.pl
if [ ${isSingleEnd} = 'TRUE' ]; then
    bam2frag=simply_bam2frags_singleEnd.pl
fi


${PERL_PATH}/perl ${curr_dir}/src/${bam2frag} --read_file ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.MAPQ${MAPQ}.bam \
        --output_file ${qc_dir}/${OUTPUT_PREFIX}.fragments.tsv --samtools_path $SAMTOOLS_PATH

${TABIX_PATH}/bgzip -f ${qc_dir}/${OUTPUT_PREFIX}.fragments.tsv.len 

## sort and index fragment file
# faster sort by R data.table
${R_PATH}/Rscript --vanilla ${curr_dir}/src/sort_frags.R ${qc_dir}/${OUTPUT_PREFIX}.fragments.tsv

## index fragment file
#sort -k1,1 -k2,2n -T ${mapRes_dir}/tmp/  ${qc_dir}/${OUTPUT_PREFIX}.fragments.tsv > ${qc_dir}/${OUTPUT_PREFIX}.fragments.sorted.tsv
#mv ${qc_dir}/${OUTPUT_PREFIX}.fragments.sorted.tsv ${qc_dir}/${OUTPUT_PREFIX}.fragments.tsv
${TABIX_PATH}/bgzip -f ${qc_dir}/${OUTPUT_PREFIX}.fragments.tsv 
${TABIX_PATH}/tabix -f -p bed ${qc_dir}/${OUTPUT_PREFIX}.fragments.tsv.gz


