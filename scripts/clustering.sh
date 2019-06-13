#!/bin/bash


mtx_file=$1

output_dir=$2 
output_dir=${output_dir}/downstream_analysis/${CELL_CALLER}
mkdir -p $output_dir

curr_dir=`dirname $0`

${R_PATH}/Rscript --vanilla ${curr_dir}/Utils/clustering.R $mtx_file $CLUSTERING_METHOD $K_CLUSTERS $output_dir $GENOME_NAME 
