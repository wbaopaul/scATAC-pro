#!/bin/bash

set -e

input_bam=$1  ## two input bam file for two cluster, seperated by ,

bams=(${input_bam//,/ })

bam1=${bams[0]}
bam2=${bams[1]}

output_dir=$2
output_dir=${output_dir}/downstream_analysis/${CELL_CALLER}/data_by_cluster/footprint
mkdir -p $output_dir

input_peak=${2}/peaks/${PEAK_CALLER}/${OUTPUT_PREFIX}.${MAPPING_METHOD}_peaks_BlacklistRemoved.narrowPeak  ## suppose path for peak file


unset PYTHONPATH

echo "predict footprint..."
cluster1=${bam1/.bam/}
cluster1=$(basename $cluster1)
${HINT_PATH}/rgt-hint footprinting --atac-seq --paired-end --organism=${GENOME_NAME} --output-location=$output_dir  \
        --output-prefix=${cluster1} $bam1 $input_peak

cluster2=${bam2/.bam/}
cluster2=$(basename $cluster2)
${HINT_PATH}/rgt-hint footprinting --atac-seq --paired-end --organism=${GENOME_NAME} --output-location=$output_dir  \
        --output-prefix=${cluster2} $bam2 $input_peak

echo "overlap with motif annotation ... "
${HINT_PATH}/rgt-motifanalysis matching --organism=${GENOME_NAME} --input-files ${output_dir}/${cluster1}.bed ${output_dir}/${cluster2}.bed --output-location $output_dir

ncore=$(nproc --all)
ncore=$(($ncore - 1))
echo "number of cores: $ncore"
echo "Differential binding analysis ..."
${HINT_PATH}/rgt-hint differential --organism=${GENOME_NAME} --bc --nc 10 --mpbs-file1=${output_dir}/${cluster1}_mpbs.bed \
    --mpbs-file2=${output_dir}/${cluster1}_mpbs.bed --reads-file1=$bam1 --reads-file2=$bam2 --condition1=$cluster1 --condition2=$cluster2 --output-location=${output_dir}/${cluster1}_${cluster2}


echo "Done!"


