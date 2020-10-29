#!/bin/bash

## reconstruct peak-cell matrix given a union peak file, a fragment.txt file and
## a barcodes.txt file, separated by comma
## the output peak-by-cell matrix will be saved under reConstruct_Matrix/, which
## was under the same directory as barcodes.txt file

set -e

inputs=$1 

inputs=(${inputs//,/ })
seurat_file=${inputs[0]}
qc_file=${inputs[1]}

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"


${R_PATH}/R --vanilla --args $seurat_file $qc_file < ${curr_dir}/src/addQC2seurat.R 

