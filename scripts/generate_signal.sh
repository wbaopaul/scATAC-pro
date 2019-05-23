#!/bin/bash

input_bam=$1

signal_dir=${2}/signal

mkdir -p $signal_dir

bamName=${input_bam##*/}   ## use input bam filename as prefix

# index bam file
${SAMTOOLS_PATH}/samtools index $input_bam

# generate bw file
unset PYTHONPATH
${DEEPTOOLS_PATH}/bamCoverage --numberOfProcessors max \
  --bam $input_bam --binSize 200 --skipNonCoveredRegions \
  --outFileName ${signal_dir}/${bamName}.bw

# generate bedgraph
${DEEPTOOLS_PATH}/bamCoverage --numberOfProcessors max \
  --bam $input_bam --binSize 200 --skipNonCoveredRegions \
  --outFileFormat bedgraph --outFileName ${signal_dir}/${bamName}.bedgraph
