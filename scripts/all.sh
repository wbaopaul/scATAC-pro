#!/bin/bash

set -e
unset PYTHONPATH
curr_dir=`dirname $0`
${curr_dir}/process.sh $1 $2 $3

#downstream
mat_file=${OUTPUT_DIR}/filtered_matrix/${CELL_CALLER}/matrix.mtx
${curr_dir}/downstream.sh $mat_file $2 $3

