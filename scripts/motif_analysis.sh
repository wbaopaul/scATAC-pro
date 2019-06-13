#!/bin/bash


mtx_file=$1

output_dir=$2 
output_dir=${output_dir}/downstream_analysis
mkdir -p $output_dir

curr_dir=`dirname $0`

${R_PATH}/Rscript --vanilla ${curr_dir}/Utils/motif_analysis.R $mtx_file $GENOME_NAME $output_dir 
