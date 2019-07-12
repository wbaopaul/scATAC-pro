#!/bin/bash

set -e
unset PYTHONPATH
curr_dir=`dirname $0`
#preprocess
#make INPUT_FILE=$1 CONFIG_FILE=$2 CONFIG_SYS=$3 preprocess > ${logDir}/preprocess.log 2>&1  
${curr_dir}/preprocess.sh $1 $2 $3

#downstream
if [ $USE_BIN ]; then
    mat_file=${OUTPUT_DIR}/raw_matrix/bin_mat_resl${BIN_RESL}/${outfile_prefix}.bin_resl${BIN_RESL}.mtx
else
    mat_file=${OUTPUT_DIR}/raw_matrix/peak_mat/${outfile_prefix}.peak.barcode.mtx
fi
#make INPUT_FILE=$mat_file CONFIG_FILE=$2 CONFIG_SYS=$3 downstream > ${logDir}/downstream.log 2>&1  
${curr_dir}/downstream.sh $mat_file $2 $3

#report
#make INPUT_FILE=$OUTPUT_DIR CONFIG_FILE=$2 CONFIG_SYS=$3 report > ${logDir}/report.log 2>&1  
${curr_dir}/report.sh ${OUTPUT_DIR}/summary $2 $3
