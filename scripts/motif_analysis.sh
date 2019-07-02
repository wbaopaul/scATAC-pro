#!/bin/bash


mtx_file=$1
source $2
output_dir=${OUTPUT_DIR}/downstream_analysis/${CELL_CALLER}
mkdir -p $output_dir

curr_dir=`dirname $0`

${R_PATH}/Rscript --vanilla ${curr_dir}/src/motif_analysis.R $mtx_file $GENOME_NAME $output_dir 
