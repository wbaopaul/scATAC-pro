#!/bin/bash

set -e

ff=$1  ## a bed or txt file, with Barcode and Cluster information

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

output_dir=${OUTPUT_DIR}/downstream_analysis/${CELL_CALLER}/data_by_cluster
mkdir -p $output_dir

input_bam=${OUTPUT_DIR}/mapping_result/${OUTPUT_PREFIX}.${MAPPING_METHOD}.positionsort.MAPQ${MAPQ}.bam  ## suppose path for bam file

curr_dir=`dirname $0`

## save header
${SAMTOOLS_PATH}/samtools view -H $input_bam > header.txt



## split bam
${PERL_PATH}/perl ${curr_dir}/src/split_bam2clusters.pl --cluster_file $ff --bam_file $input_bam \
    --output_dir $output_dir --samtools_path $SAMTOOLS_PATH

unset PYTHONPATH
ncore=$(nproc --all)
ncore=$(($ncore - 1))
for file0 in $(find $output_dir -name *cluster*.sam); do
    fname_bam=${file0/.sam/.bam}
    cat header.txt $file0 > ${file0}.wheader
    ${SAMTOOLS_PATH}/samtools view -bS -@ $ncore ${file0}.wheader > $fname_bam 
    rm ${file0}.wheader
    rm $file0
    echo "generate bw file..."
    fname_bw=${file0/.sam/.bw}
    ${SAMTOOLS_PATH}/samtools index -@ $ncore $fname_bam
    ${DEEPTOOLS_PATH}/bamCoverage --numberOfProcessors max \
         --bam $fname_bam --binSize 10 --skipNonCoveredRegions \
         --outFileName $fname_bw
done





