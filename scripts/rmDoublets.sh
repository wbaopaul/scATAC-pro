#!/bin/bash


inputs=$1  ## seurat obj,and expeted doublet rate,
 
# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"


inputs=(${inputs//,/ })
seuratObj_file=${inputs[0]}
drate=${inputs[1]}

output_dir=`dirname $seuratObj_file`
curr_dir=`dirname $0`

${R_PATH}/Rscript --vanilla ${curr_dir}/src/rmDoublets.R $seuratObj_file $drate

echo "Doublets removed!"
