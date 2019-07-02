#!/bin/bash


mtx_file=$1
source $2
 
output_dir=${OUTPUT_DIR}/downstream_analysis/${CELL_CALLER}
mkdir -p $output_dir

curr_dir=`dirname $0`

${R_PATH}/Rscript --vanilla ${curr_dir}/src/clustering.R $mtx_file $CLUSTERING_METHOD $K_CLUSTERS $output_dir $GENOME_NAME 
