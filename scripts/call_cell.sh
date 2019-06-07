#!/bin/bash

### will search peak file in peaks/*narrowPeak under $2
set -e

input_mtx_file=$1
ABS_PATH=`cd "$2"; pwd`

bc_stat_file=${ABS_PATH}/qc_result/${OUTPUT_PREFIX}.${MAPPING_METHOD}.qc_per_barcode.bed
mapping_qc_file=${ABS_PATH}/qc_result/${OUTPUT_PREFIX}.${MAPPING_METHOD}.summary
fragments_file=${ABS_PATH}/raw_matrix/fragments.bed
output_dir=${ABS_PATH}/filtered_matrix

curr_dir=`dirname $0`


echo "${CELL_CALLER} is usded for cell calling..."


if [ ${CELL_CALLER} = 'EmptyDrop' ];then
	output_dir1=${output_dir}/EmptyDrop
	mkdir -p $output_dir1
	${R_PATH}/R --vanilla --args $input_mtx_file $output_dir1 ${EmptyDrop_FDR} < ${curr_dir}/EmptyDrop.R
	
	echo "Now reporting barcode qc stat "
        ABS_PATH=`cd "$output_dir1"; pwd`
  
        ${R_PATH}/R --vanilla --args ${ABS_PATH}/report_kneepoint.html $mapping_qc_file $bc_stat_file \
        ${ABS_PATH}/kneepoint/barcodes.txt $fragments_file $curr_dir < ${curr_dir}/render2report.R
        
	${R_PATH}/R --vanilla --args ${ABS_PATH}/report_fdr${EmptyDrop_FDR}.html $mapping_qc_file $bc_stat_file \
        ${ABS_PATH}/fdr${EmptyDrop_FDR}/barcodes.txt $fragments_file $curr_dir  < ${curr_dir}/render2report.R

elif [ ${CELL_CALLER} = 'FILTER' ];then
	output_dir1=${output_dir}/FILTER
	mkdir -p $output_dir1
        ${R_PATH}/Rscript  ${curr_dir}/filter_barcodes.R  --bc_stat_file $bc_stat_file --raw_mtx_file $input_mtx_file \
        --output_prefix $output_dir1 $BC_FILTER_CUTOFF
        
	echo "Now reporting barcode qc stat "
    	ABS_PATH=`cd "$output_dir1"; pwd`
    
    	${R_PATH}/R --vanilla --args ${ABS_PATH}/report.html $mapping_qc_file $bc_stat_file \
        ${ABS_PATH}/barcodes.txt $fragments_file $curr_dir < ${curr_dir}/render2report.R

else 
	## using cellranger
	## calculate peak coverage fraction
        peak_file=${2}/peaks/MACS2/${OUTPUT_PREFIX}.${MAPPING_METHOD}_peaks_BlacklistRemoved.narrowPeak
	gsize=$(awk '{sum+=$2}; END {print sum}' $CHROM_SIZE_FILE)
	peaks_esize=$(awk '{sum+=$3}; END {print sum}' $peak_file)
	peaks_ssize=$(awk '{sum+=$2}; END {print sum}' $peak_file)
	peak_cov=$(( ${peaks_esize} - ${peaks_ssize} ))
	peak_cov_frac=$(${PERL_PATH}/perl -E "say $peak_cov/$gsize")

	## implement cellrange cell calling
	output_dir1=${output_dir}/cellranger
	mkdir -p output_dir1
	
        ${R_PATH}/R --vanilla --args $input_mtx_file $output_dir1 $peak_cov_frac $bc_stat_file \
        < ${curr_dir}/cellranger_cell_caller.R
    
    	echo "Now reporting barcode qc stat "
    	ABS_PATH=`cd "$output_dir1"; pwd`
    
    	${R_PATH}/R --vanilla --args ${ABS_PATH}/report.html $mapping_qc_file $bc_stat_file \
	${ABS_PATH}/barcodes.txt $fragments_file $curr_dir < ${curr_dir}/render2report.R

fi



