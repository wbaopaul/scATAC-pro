#!/bin/bash

set -e
unset PYTHONPATH

input_mtx=$1

curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf $2
read_conf $3

outfile_prefix=${OUTPUT_PREFIX}.${MAPPING_METHOD}

frag_file=${OUTPUT_DIR}/summary/fragments.bed

## clustering
${curr_dir}/clustering.sh $input_mtx $2 $3 &

## motif analysis
${curr_dir}/motif_analysis.sh $input_mtx $2 $3 &

wait

## do DA
${curr_dir}/doDA.sh $input_mtx $2 $3 &

SPLIT_BAM2CLUSTER=${SPLIT_BAM2CLUSTER^^}
## split bam to cluster
if [ "$SPLIT_BAM2CLUSTER" = "TRUE" ]; then
    input_cluster_table=${OUTPUT_DIR}/downstream_analysis/${CELL_CALLER}/cell_cluster_table.txt 
    ${curr_dir}/split_bam2clusters.sh $input_cluster_table $2 $3 &
fi

wait

## footprinting analysis
DO_FOOTPRINT=${DO_FOOTPRINT^^}
if [ "$DO_FOOTPRINT" = "TRUE" ]; then
    bam1=${OUTPUT_DIR}/downstream_analysis/${CELL_CALLER}/data_by_cluster/cluster_${cluster1}.bam
    bam2=${OUTPUT_DIR}/downstream_analysis/${CELL_CALLER}/data_by_cluster/cluster_${cluster2}.bam
    ${curr_dir}/footprint.sh ${bam1},${bam2} $2 $3 
fi

${curr_dir}/report.sh ${OUTPUT_DIR}/summary $2 $3

