#!/bin/bash

set -e
unset PYTHONPATH

input_mtx=$1

curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf $2
read_conf $3

outfile_prefix=${OUTPUT_PREFIX}.${MAPPING_METHOD}

frag_file=${OUTPUT_DIR}/raw_matrix/fragments.bed

## clustering
make INPUT_FILE=$input_mtx CONFIG_FILE=$2 CONFIG_SYS=$3 clustering > ${logDir}/clustering.log 2>&1 & 

## motif analysis
make INPUT_FILE=$input_mtx CONFIG_FILE=$2 CONFIG_SYS=$3 motif_analysis > ${logDir}/motif_analysis.log 2>&1 & 

## split bam to cluster
if [ $SPLIT_BAM2CLUSTER ]; then
    input_cluster_table=${OUTPUT_DIR}/downstream_analysis/${CELL_CALLER}/cell_cluster_table.txt 
    make INPUT_FILE=$input_cluster_table CONFIG_FILE=$2 CONFIG_SYS=$3 split_bam > ${logDir}/split_bam.log 2>&1 & 
fi

## footprinting analysis
if [ $DO_FOOTPRINT ]; then
    bam1=${OUTPUT_DIR}/downstream_analysis/${CELL_CALLER}/data_by_cluster/${cluster1}.bam
    bam2=${OUTPUT_DIR}/downstream_analysis/${CELL_CALLER}/data_by_cluster/${cluster2}.bam
    make INPUT_FILE=${bam1},${bam2} CONFIG_FILE=$2 CONFIG_SYS=$3 footprint > ${logDir}/footprint.log 2>&1 & 
fi



