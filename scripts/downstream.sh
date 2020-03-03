#!/bin/bash

set -e
unset PYTHONPATH

input_mtx=$1

curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf $2
read_conf $3


frag_file=${OUTPUT_DIR}/summary/${OUTPUT_PREFIX}.fragments.txt

## clustering
${curr_dir}/clustering.sh $input_mtx $2 $3 

## motif analysis
${curr_dir}/motif_analysis.sh $input_mtx $2 $3 

seurat_obj=${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}/seurat_obj.rds
## do DA
if [ "$RUN_DA" = "TRUE" ]; then
    ${curr_dir}/runDA.sh $seurat_obj $2 $3 &
fi

SPLIT_BAM2CLUSTER=${SPLIT_BAM2CLUSTER^^}
## split bam to cluster
if [ "$SPLIT_BAM2CLUSTER" = "TRUE" ]; then
    input_cluster_table=${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}/cell_cluster_table.txt 
    ${curr_dir}/split_bam2clusters.sh $input_cluster_table $2 $3 &
fi

wait

## go analysis
if [ "$RUN_GO" = "TRUE" ]; then
    ${curr_dir}/runGO.sh ${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}/differential_accessible_features_${group1}_vs_${group2}.txt $2 $3 &
fi

if [ "$RUN_Cicero" = "TRUE" ]; then
    ${curr_dir}/runCicero.sh $seurat_obj $2 $3 &
fi

## footprinting analysis
${curr_dir}/footprint.sh ${group1_fp},${group2_fp} $2 $3

${curr_dir}/report.sh ${OUTPUT_DIR}/summary $2 $3

