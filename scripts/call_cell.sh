#!/bin/bash

### will search peak file in peaks/*narrowPeak under $2
set -e

input_mtx_file=$1
output_dir=${2}/filtered_matrix

if [ $cell_caller == 'EmptyDrop' ];then
	output_prefix=${output_dir}/EmptyDrop_mat
	${R_PATH}/R --vanilla --args $input_mtx_file $output_prefix ${EmptyDrop_fdr} < EmptyDrop.R
fi
