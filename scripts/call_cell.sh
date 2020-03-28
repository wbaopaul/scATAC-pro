#!/bin/bash

### will search peak file in peaks/*narrowPeak under $2
set -e

input_mtx_file=$1
# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"
ABS_PATH=`cd "$OUTPUT_DIR"; pwd`

bc_stat_file=${ABS_PATH}/summary/${OUTPUT_PREFIX}.${PEAK_CALLER}.qc_per_barcode.txt
mapping_qc_file=${ABS_PATH}/summary/${OUTPUT_PREFIX}.MappingStats
fragments_file=${ABS_PATH}/summary/${OUTPUT_PREFIX}.fragments.txt
output_dir=${ABS_PATH}/filtered_matrix

curr_dir=`dirname $0`


echo "${CELL_CALLER} is used for cell calling..."


if [ ${CELL_CALLER} = 'EmptyDrop' ];then
	output_dir1=${output_dir}/${PEAK_CALLER}/EmptyDrop
	mkdir -p $output_dir1
	${R_PATH}/R --vanilla --args $input_mtx_file $output_dir1 ${EmptyDrop_FDR} < ${curr_dir}/src/EmptyDrop.R
	
  

elif [ ${CELL_CALLER} = 'FILTER' ];then
	output_dir1=${output_dir}/${PEAK_CALLER}/FILTER
	mkdir -p $output_dir1
    
    echo "${R_PATH}/Rscript  ${curr_dir}/src/filter_barcodes.R  --bc_stat_file $bc_stat_file --raw_mtx_file $input_mtx_file --output_dir $output_dir1 ${FILTER_BC_CUTOFF}" > ${output_dir1}/tmpJob
    bash ${output_dir1}/tmpJob
    rm ${output_dir1}/tmpJob
        
    

else 
	## using cellranger
	## calculate peak coverage fraction
	gsize=$(awk '{sum+=$2}; END {print sum}' $CHROM_SIZE_FILE)

	## implement cellrange cell calling
	output_dir1=${output_dir}/${PEAK_CALLER}/cellranger
	mkdir -p $output_dir1
	
    ${R_PATH}/R --vanilla --args $input_mtx_file $output_dir1 $gsize $bc_stat_file \
    < ${curr_dir}/src/cellranger_cell_caller.R

fi


## generate report
#${curr_dir}/report.sh ${OUTPUT_DIR}/summary $2 $3


