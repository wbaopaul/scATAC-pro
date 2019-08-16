#!/bin/bash

input_bam=$1
# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

signal_dir=${OUTPUT_DIR}/signal

mkdir -p $signal_dir


# index bam file
ncore=`nproc --all`
ncore=$(($ncore - 1))

if [[ ! -f ${input_bam}.bai ]];then
    ${SAMTOOLS_PATH}/samtools index -@ $ncore ${input_bam}.bam
fi

echo "generate bw file..."
unset PYTHONPATH
${DEEPTOOLS_PATH}/bamCoverage --numberOfProcessors max \
  --bam $input_bam --binSize 20 --skipNonCoveredRegions --normalizeUsing BPM \
  --outFileName ${signal_dir}/${OUTPUT_PREFIX}.aggregated.bw &

#echo "generate bedgraph..."
${DEEPTOOLS_PATH}/bamCoverage --numberOfProcessors max \
  --bam $input_bam --binSize 20 --skipNonCoveredRegions -- --normalizeUsing BPM \
  --outFileFormat bedgraph --outFileName ${signal_dir}/${OUTPUT_PREFIX}.aggregated.bedgraph &
wait

echo "generate count around TSS..."
${DEEPTOOLS_PATH}/computeMatrix reference-point -S ${signal_dir}/${OUTPUT_PREFIX}.aggregated.bw -R $TSS \
    -a 1200 -b 1200 -o ${signal_dir}/${OUTPUT_PREFIX}.aggregated.mtx.gz



