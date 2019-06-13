#!/bin/bash


report_dir=$1
mkdir -p $report_dir
abs_report_dir=`cd ${report_dir}; pwd`
echo "report path: ${abs_report_dir}"

abs_out_dir=`cd ${2}; pwd`

downstream_dir=${abs_out_dir}/downstream_analysis/${CELL_CALLER}
mkdir -p $downstream_dir

mapping_qc_file=${abs_out_dir}/qc_result/${OUTPUT_PREFIX}.${MAPPING_METHOD}.summary
bc_stat_file=${abs_out_dir}/qc_result/${OUTPUT_PREFIX}.${MAPPING_METHOD}.qc_per_barcode.bed
barcode_file=${abs_out_dir}/filtered_matrix/${CELL_CALLER}/barcodes.txt
fragment_file=${abs_out_dir}/fragments/fragments_MAPQ${MAPQ}.bed
tss_escore_file=${abs_out_dir}/signal/${OUTPUT_PREFIX}.${MAPPING_METHOD}.dedup.paired.MAPQ${MAPQ}.bam.mtx

curr_dir=`dirname $0`


${R_PATH}/Rscript --vanilla ${curr_dir}/Utils/render2report.R \
    ${abs_report_dir}/scATAC-pro_report_${CELL_CALLER}.html $mapping_qc_file $bc_stat_file \
    $barcode_file $fragments_file $downstream_dir $tss_escore_file
         

