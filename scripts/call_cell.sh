#!/bin/bash

### will search peak file in peaks/*narrowPeak under $2
set -e

input_mtx_file=$1
ABS_PATH=`cd "$2"; pwd`

bc_stat_file=${ABS_PATH}/qc_result/${OUTPUT_PREFIX}.${MAPPING_METHOD}.qc_per_barcode.bed
output_dir=${ABS_PATH}/filtered_matrix
curr_dir=`dirname $0`

echo "${CELL_CALLER} is usded for cell calling..."


if [ ${CELL_CALLER} = 'EmptyDrop' ];then
	output_prefix=${output_dir}/EmptyDrop_mat/
	mkdir -p $output_prefix
	${R_PATH}/R --vanilla --args $input_mtx_file $output_prefix ${EmptyDrop_FDR} < ${curr_dir}/EmptyDrop.R
	
	echo "Now reporting barcode qc stat "
    ABS_PATH=`cd "$output_prefix"; pwd`
  
    ${R_PATH}/R --vanilla --args ${ABS_PATH}/bc_stat_kneepoint.html ${bc_stat_file} ${ABS_PATH}/filtered.kneepoint.barcodes.txt $curr_dir \
    < ${curr_dir}/render2report.R
     ${R_PATH}/R --vanilla --args ${ABS_PATH}/bc_stat_fdr${EmptyDrop_FDR}.html ${bc_stat_file} ${ABS_PATH}/filtered.${EmptyDrop_FDR}.barcodes.txt \          $curr_dir  < ${curr_dir}/render2report.R

elif [ ${CELL_CALLER} = 'FILTER' ];then
	output_prefix=${output_dir}/FILTER_mat/
	mkdir -p $output_prefix
    ${R_PATH}/Rscript  ${curr_dir}/filter_barcodes.R  --bc_stat_file $bc_stat_file --raw_mtx_file $input_mtx_file \
    --output_prefix $output_prefix $BC_CUTOFF
    echo "Now reporting barcode qc stat "
    ABS_PATH=`cd "$output_prefix"; pwd`
    
    ${R_PATH}/R --vanilla --args ${ABS_PATH}/bc_stat.html ${bc_stat_file} ${ABS_PATH}/filtered.barcodes.txt $curr_dir < ${curr_dir}/render2report.R

else 
	## using cellranger
	## calculate peak coverage fraction
    peak_file=${2}/peaks/${OUTPUT_PREFIX}.${MAPPING_METHOD}_peaks_filterBlacklist.narrowPeak
	gsize=$(awk '{sum+=$2}; END {print sum}' $CHROM_SIZE_FILE)
	peaks_esize=$(awk '{sum+=$3}; END {print sum}' $peak_file)
	peaks_ssize=$(awk '{sum+=$2}; END {print sum}' $peak_file)
	peak_cov=$(( ${peaks_esize} - ${peaks_ssize} ))
	peak_cov_frac=$(${PERL_PATH}/perl -E "say $peak_cov/$gsize")

	## implement cellrange cell calling
	output_prefix=${output_dir}/cellranger_mat/
	mkdir -p output_prefix
	
    ${R_PATH}/R --vanilla --args $input_mtx_file $output_prefix $peak_cov_frac $bc_stat_file  < ${curr_dir}/cellranger_cell_caller.R
    
    echo "Now reporting barcode qc stat "
    ABS_PATH=`cd "$output_prefix"; pwd`
    
    ${R_PATH}/R --vanilla --args ${ABS_PATH}/bc_stat.html ${bc_stat_file} ${ABS_PATH}/filtered.barcodes.txt $curr_dir < ${curr_dir}/render2report.R

fi



