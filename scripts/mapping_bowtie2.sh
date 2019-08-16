#!/bin/bash



input_ff=$1
fastqs=(${input_ff//,/ })
# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

mapRes_dir="${OUTPUT_DIR}/mapping_result"


echo "Starting bowtie2 alignment ... "
${BOWTIE2_PATH}/bowtie2 -1 ${fastqs[0]} -2 ${fastqs[1]} -x $BOWTIE2_INDEX  $BOWTIE2_OPTS -S ${mapRes_dir}/${OUTPUT_PREFIX}.sam 

echo "Bowtie2 Mapping Done!"


## convert to bam
echo "Converting sam to bam ... "
ncore=$(nproc --all)
ncore=$(($ncore - 1))

${SAMTOOLS_PATH}/samtools view -@ $ncore -h -bS ${mapRes_dir}/${OUTPUT_PREFIX}.sam > ${mapRes_dir}/${OUTPUT_PREFIX}.bam

rm ${mapRes_dir}/${OUTPUT_PREFIX}.sam

