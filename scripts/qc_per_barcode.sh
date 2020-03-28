#!/bin/bash

set -e
inputs=$1
inputs=(${inputs//,/ })
input_frags=${inputs[0]}  ## fragment.txt file
peaks=${inputs[1]}  ## peak file

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

qc_dir="${OUTPUT_DIR}/summary"
mkdir -p $qc_dir



## use perl script to get the matrix
curr_dir=`dirname $0`

file_prefix=${OUTPUT_PREFIX}

#peaks=${OUTPUT_DIR}/peaks/${PEAK_CALLER}/${file_prefix}_features_BlacklistRemoved.bed


${R_PATH}/R --vanilla --args $input_frags $peaks $PROMOTERS $TSS $ENHANCERS ${qc_dir}/${file_prefix}.${PEAK_CALLER}.qc_per_barcode.txt < ${curr_dir}/src/get_qc_per_barcode.R


echo "QC per barcode done! "



