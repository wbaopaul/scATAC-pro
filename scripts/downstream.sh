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
${curr_dir}/clustering.sh $input_mtx $2 $3 &

## motif analysis
${curr_dir}/motif_analysis.sh $input_mtx $2 $3 &

wait
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
    ${curr_dir}/runGO.sh ${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}/differential_accessible_features_cluster_${group1}_VS_cluster_${group2}.txt $2 $3 &
fi

if [ "$RUN_Cicero" = "TRUE" ]; then
    ${curr_dir}/runCicero.sh $seurat_obj $2 $3 &
fi

## footprinting analysis
bam_file=${OUTPUT_DIR}/mapping_result/cell_barcodes.positionsort.MAPQ${MAPQ}.bam
cls=($(cut -f2 $input_cluster_table | sed '1d'| sort -u))
DO_FOOTPRINT=${DO_FOOTPRINT^^}
if [ "$DO_FOOTPRINT" = "TRUE" ]; then
    groups1=(${group1_fp/:/ })
    groups2=(${group2_fp/:/ })

    if [ "$group1_fp" == 'one' ]; then
        echo "do all one-vs-rest comparison, this would take a very long time..."
        for cl0 in "${cls[@]}"
        do
            bam1=${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}/data_by_cluster/cluster_${cl0}.bam
            ${curr_dir}/footprint_by_cluster.sh ${bam1},${bam_file} $2 $3 
        done
    elif [ "$group2_fp" == 'rest' ]; then
        for cl0 in "${groups1[@]}"
        do
            bam1=${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}/data_by_cluster/cluster_${cl0}.bam
            ${curr_dir}/footprint_by_cluster.sh ${bam1},${bam_file} $2 $3 
        done
    else
        bam1=${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}/data_by_cluster/cluster_${group1_fp}.bam
        bam2=${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}/data_by_cluster/cluster_${group2_fp}.bam
        ${curr_dir}/footprint_by_cluster.sh ${bam1},${bam2} $2 $3 
    fi
fi


${curr_dir}/report.sh ${OUTPUT_DIR}/summary $2 $3

