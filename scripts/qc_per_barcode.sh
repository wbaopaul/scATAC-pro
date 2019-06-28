#!/bin/bash

set -e

input_frags=$1  ## use read pair(fragment) bed file as input
source $2
qc_dir="${OUTPUT_DIR}/qc_result"
mkdir -p $qc_dir



## use perl script to get the matrix
curr_dir=`dirname $0`

file_prefix=${OUTPUT_PREFIX}.${MAPPING_METHOD}

peaks=${2}/peaks/${PEAK_CALLER}/${file_prefix}_peaks_BlacklistRemoved.narrowPeak


${R_PATH}/R --vanilla --args $input_frags $peaks $PROMOTERS $TSS $ENHANCERS ${qc_dir}/${file_prefix}.qc_per_barcode.bed < ${curr_dir}/Utils/get_qc_per_barcode.R


echo "QC per barcode done! "



