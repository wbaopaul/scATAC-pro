#!/bin/bash

set -e

input_frags=$1  ## use read pair(fragment) bed file as input

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

qc_dir="${OUTPUT_DIR}/summary"
mkdir -p $qc_dir



## use perl script to get the matrix
curr_dir=`dirname $0`

file_prefix=${OUTPUT_PREFIX}.${MAPPING_METHOD}

peaks=${OUTPUT_DIR}/peaks/${PEAK_CALLER}/${file_prefix}_features_BlacklistRemoved.bed


${R_PATH}/R --vanilla --args $input_frags $peaks $PROMOTERS $TSS $ENHANCERS ${qc_dir}/${file_prefix}.qc_per_barcode.bed < ${curr_dir}/src/get_qc_per_barcode.R


echo "QC per barcode done! "



