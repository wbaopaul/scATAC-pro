#!/bin/bash

set -e

ff=$1  ## a bed or txt file, with Barcode and Cluster information
output_dir=$2
output_dir=${output_dir}/downstream_analysis/${CELL_CALLER}/data_by_cluster
mkdir -p $output_dir

input_bam=${2}/filtered_bam/${OUTPUT_PREFIX}.${MAPPING_METHOD}.dedup.MAPQ30.bam  ## suppose path for bam file

curr_dir=`dirname $0`

ncore=$(nproc --all)
ncore=$(($ncore - 1))
## save header
${SAMTOOLS_PATH}/samtools view -H $input_bam > header.txt

## get sam
${SAMTOOLS_PATH}/samtools view -@ $ncore $input_bam > tmp.sam

## split sam
${PERL_PATH}/perl ${curr_dir}/Utils/split_sam2clusters.pl --cluster_file $ff --sam_file tmp.sam --output_dir $output_dir

rm tmp.sam

unset PYTHONPATH
## cat head to each sam and tranfer back to bam
for file0 in $(find $output_dir -name *cluster*.sam); do
    cat header.txt $file0 > ${file0}_withheader
    fname_bam=${file0/.sam/.bam}
    ${SAMTOOLS_PATH}/samtools view -@ $ncore -bS -h ${file0}_withheader > $fname_bam
    rm ${file0}_withheader
    rm $file0
    echo "generate bw file..."
    fname_bw=${fname_bam/.bam/.bw}
    ${SAMTOOLS_PATH}/samtools index -@ $ncore $fname_bam
    ${DEEPTOOLS_PATH}/bamCoverage --numberOfProcessors max \
         --bam $fname_bam --binSize 1 --skipNonCoveredRegions \
         --outFileName $fname_bw
done





