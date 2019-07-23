#!/bin/bash



input_ff=$1
fastqs=(${input_ff//,/ })

mapRes_dir="${2}"


echo "Starting bowtie2 alignment ... "
${BOWTIE2_PATH}/bowtie2 -1 ${fastqs[0]} -2 ${fastqs[1]} -x $BOWTIE2_INDEX  $BOWTIE2_OPTS -S ${mapRes_dir}/${OUTPUT_PREFIX}.bowtie2.sam 

echo "Bowtie2 Mapping Done!"


## convert to bam
echo "Converting sam to bam ... "
ncore=$(nproc --all)
ncore=$(($ncore - 1))

${SAMTOOLS_PATH}/samtools view -@ $ncore -h -bS ${mapRes_dir}/${OUTPUT_PREFIX}.bowtie2.sam > ${mapRes_dir}/${OUTPUT_PREFIX}.bowtie2.bam

rm ${mapRes_dir}/${OUTPUT_PREFIX}.bowtie2.sam

