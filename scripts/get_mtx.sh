#!/bin/bash

### will search fragments.bed under output_dir/summary
set -e
feature_file=$1  ## peak file

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

mat_dir="${OUTPUT_DIR}/raw_matrix"
mkdir -p $mat_dir

curr_dir=`dirname $0`


echo "Getting peak by barcode matrix..."
    # this R script will output sorted fragments as well
    mtx_peak_dir=${mat_dir}/${PEAK_CALLER}
    mkdir -p ${mtx_peak_dir}
    ${R_PATH}/R --vanilla --args ${OUTPUT_DIR}/summary/fragments.bed $feature_file ${mtx_peak_dir} 5 5 < ${curr_dir}/src/get_mtx.R 


echo "Get Matrix Done!"






