#!/bin/bash


mtx_file=$1

output_dir=$2 
output_dir=${output_dir}/downstream_analysis/${CELL_CALLER}
mkdir -p $output_dir

curr_dir=`dirname $0`


${R_PATH}/Rscript ${curr_dir}/Utils/render2report.R ${ABS_PATH}/scATAC-pro_report.html $mapping_qc_file $bc_stat_file \
             ${ABS_PATH}/barcodes.txt $fragments_file $curr_dir 
         

