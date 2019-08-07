#!/bin/bash

set -e
unset PYTHONPATH
curr_dir=`dirname $0`
#preprocess
#make INPUT_FILE=$1 CONFIG_FILE=$2 CONFIG_SYS=$3 preprocess > ${logDir}/preprocess.log 2>&1  
${curr_dir}/preprocess.sh $1 $2 $3

#downstream
    mat_file=${OUTPUT_DIR}/raw_matrix/${PEAK_CALLER}/matrix.mtx
${curr_dir}/downstream.sh $mat_file $2 $3

