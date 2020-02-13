#!/bin/bash

set -e

ff=$1  ## a bed or txt file, with Barcode and Cluster information

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

output_dir=${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}/data_by_cluster
mkdir -p $output_dir

input_bam=${OUTPUT_DIR}/mapping_result/${OUTPUT_PREFIX}.positionsort.MAPQ${MAPQ}.bam  ## suppose path for bam file

curr_dir=`dirname $0`



## split bam
${PERL_PATH}/perl ${curr_dir}/src/split_bam2clusters.pl --cluster_file $ff --bam_file $input_bam \
    --output_dir $output_dir --samtools_path $SAMTOOLS_PATH

unset PYTHONPATH
ncore=$(nproc --all)
ncore=$((${ncore}/2))
for file0 in $(find $output_dir -name "*cluster*.bam")
do
    echo "generate bw file..."
    fname_bw=${file0/.bam/.bw}
    ${SAMTOOLS_PATH}/samtools index -@ $ncore $file0
    ${DEEPTOOLS_PATH}/bamCoverage --numberOfProcessors max --normalizeUsing BPM \
         --bam $file0 --binSize 100 --skipNonCoveredRegions \
         --outFileName $fname_bw  &
    fname_bedgraph=${file0/.bam/.bedgraph}
    ${DEEPTOOLS_PATH}/bamCoverage --numberOfProcessors max --normalizeUsing BPM \
         --outFileFormat bedgraph --bam $file0 --binSize 100 --skipNonCoveredRegions \
         --outFileName $fname_bedgraph &
    wait
done

## move signal to output/signal
mv ${output_dir}/*bw ${OUTPUT_DIR}/signal/
mv ${output_dir}/*bedgraph ${OUTPUT_DIR}/signal/

echo "The bam file was split into different clusters!"



