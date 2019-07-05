#!/bin/bash


report_dir=$1
mkdir -p $report_dir

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

abs_report_dir=`cd ${report_dir}; pwd`
echo "report path: ${abs_report_dir}"

abs_out_dir=`cd ${OUTPUT_DIR}; pwd`

downstream_dir=${abs_out_dir}/downstream_analysis/${CELL_CALLER}
mkdir -p $downstream_dir

mapping_qc_file=${abs_out_dir}/qc_result/${OUTPUT_PREFIX}.${MAPPING_METHOD}.summary
bc_stat_file=${abs_out_dir}/qc_result/${OUTPUT_PREFIX}.${MAPPING_METHOD}.qc_per_barcode.bed
barcode_file=${abs_out_dir}/filtered_matrix/${CELL_CALLER}/barcodes.txt
fragments_file=${abs_out_dir}/raw_matrix/fragments.bed
tss_escore_file=${abs_out_dir}/signal/${OUTPUT_PREFIX}.${MAPPING_METHOD}.aggregate.mtx

curr_dir=`dirname $0`


${R_PATH}/Rscript --vanilla ${curr_dir}/src/render2report.R \
    ${abs_report_dir}/scATAC-pro_report_${CELL_CALLER}.html $mapping_qc_file $bc_stat_file \
    $barcode_file $fragments_file $downstream_dir $tss_escore_file
         

