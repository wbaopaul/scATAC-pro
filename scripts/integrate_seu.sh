#!/bin/bash


mtx_files=$1  ## mtx for each sample, separated by ,
 
# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

output_dir=${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}
mkdir -p $output_dir

curr_dir=`dirname $0`

${R_PATH}/Rscript --vanilla ${curr_dir}/src/integrate_seu.R $mtx_files $CLUSTERING_METHOD $K_CLUSTERS $output_dir $GENOME_NAME $norm_by $REDUCTION $nREDUCTION $Top_Variable_Features 
