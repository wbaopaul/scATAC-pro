#!/bin/bash

## reconstruct peak-cell matrix given a union peak file, a fragment.txt file and
## a barcodes.txt file, separated by comma
## the output peak-by-cell matrix will be saved under reConstruct_Matrix/, which
## was under the same directory as barcodes.txt file

set -e

inputs=$1 

inputs=(${inputs//,/ })
peak_file=${inputs[0]}
frag_file=${inputs[1]}
bc_file=${inputs[2]}

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

outMat_dir=`dirname $bc_file`
outMat_dir=${outMat_dir}/reConstruct_matrix
mkdir -p ${outMat_dir}

echo "re-constructing peak-by-cell matrix for each sample ..."


${R_PATH}/R --vanilla --args $frag_file $peak_file ${bc_file} ${outMat_dir} < ${curr_dir}/src/reConstMtx.R 

