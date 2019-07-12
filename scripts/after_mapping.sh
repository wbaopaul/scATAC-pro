#!/bin/bash

set -e
unset PYTHONPATH

bam_file=$1

curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf $2
read_conf $3

outfile_prefix=${OUTPUT_PREFIX}.${MAPPING_METHOD}


## 4.call peak
bam_file=${OUTPUT_DIR}/mapping_result/${outfile_prefix}.positionsort.MAPQ${MAPQ}.bam
${curr_dir}/call_peak.sh $bam_file $2 $3

## 5.generate matrix
if [ ! $USE_BIN ]; then
    feature_file=${OUTPUT_DIR}/peaks/${PEAK_CALLER}/${outfile_prefix}_peaks_BlacklistRemoved.bed
else
    feature_file="xx"
fi
${curr_dir}/get_mtx.sh $feature_file $2 $3 &

## 6.generate aggregated signal
${curr_dir}/generate_signal.sh $bam_file $2 $3 &
wait


## qc and call cell
if [ $USE_BIN ]; then
    mat_file=${OUTPUT_DIR}/raw_matrix/bin_mat_resl${BIN_RESL}/${outfile_prefix}.bin_resl${BIN_RESL}.mtx
else
    mat_file=${OUTPUT_DIR}/raw_matrix/peak_mat/${outfile_prefix}.peak.barcode.mtx
fi

frag_file=${OUTPUT_DIR}/raw_matrix/fragments.bed

if [ "$CELL_CALLER" != "filtering" ]; then
    ## 7.qc per barcode
    ${curr_dir}/qc_per_barcode.sh $frag_file $2 $3 &
    ## 8.call cell

    ${curr_dir}/call_cell.sh $mat_file $2 $3 &
else
    ## 7.qc per barcode
    
    ${curr_dir}/qc_per_barcode.sh $frag_file $2 $3 
    ## 8.call cell

    ${curr_dir}/call_cell.sh $mat_file $2 $3 
fi

wait
## report preprocessing QC
${curr_dir}/report.sh $OUTPUT_DIR/summary $2 $3
