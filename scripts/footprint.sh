#!/bin/bash

set -e

inputs=$1  ## '0,1' or one-vs-rest model like '0,rest' or 'one,rest'

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

input_cluster_table=${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}/cell_cluster_table.txt

bam_file=${OUTPUT_DIR}/mapping_result/cell_barcodes.MAPQ${MAPQ}.bam

cls=($(cut -f2 $input_cluster_table | sed '1d'| sort -u))

tmp=(${inputs//,/ })
grs1_fp=${tmp[0]}
grs2_fp=${tmp[1]}
grs1=(${grs1_fp/:/ })
grs2=(${grs2_fp/:/ })

if [ "$grs1_fp" == 'one' ]; then
    echo "do all one-vs-rest comparison, this would take a very long time..."
    for cl0 in "${cls[@]}"
    do
        echo "working on $cl0 ..."
        bam1=${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}/data_by_cluster/cluster_${cl0}.bam
        ${curr_dir}/footprint_by_cluster.sh ${bam1},${bam_file} $2 $3
    done
elif [ "$grs2_fp" == 'rest' ]; then
    for cl0 in "${grs1[@]}"
    do
        echo "working on $cl0 ..."
        bam1=${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}/data_by_cluster/cluster_${cl0}.bam
        ${curr_dir}/footprint_by_cluster.sh ${bam1},${bam_file} $2 $3
    done
else
    bam1=${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}/data_by_cluster/cluster_${grs1_fp}.bam
    bam2=${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}/data_by_cluster/cluster_${grs2_fp}.bam
    ${curr_dir}/footprint_by_cluster.sh ${bam1},${bam2} $2 $3
fi

echo "Footprint analysis done!"
