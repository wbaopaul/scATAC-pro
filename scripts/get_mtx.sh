#!/bin/bash

### will search fragments.bed under output_dir/summary
set -e
feature_file=$1  ## peak file

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

mat_dir="${OUTPUT_DIR}/raw_matrix"
mkdir -p $mat_dir


ncore=$(nproc --all)
ncore=$(($ncore - 1))



curr_dir=`dirname $0`

#echo "Sorting the output by position" -- slow in shell, just do it by R will be much faster
#sort -k1,1 -k2,2n -k3,3n ${mat_dir}/fragments.bed > ${mat_dir}/fragments.bed.sorted

USE_BIN=`echo $USE_BIN | tr a-z A-Z`

if [ "$USE_BIN" = "FALSE" ]; then 
    echo "Getting peak by barcode matrix..."
    # this R script will output sorted fragments as well
    mtx_peak_dir=${mat_dir}/peak_mat
    mkdir -p ${mtx_peak_dir}
    ${R_PATH}/R --vanilla --args ${OUTPUT_DIR}/summary/fragments.bed $feature_file ${mtx_peak_dir} < ${curr_dir}/src/get_mtx.R 
fi

if [ $USE_BIN = "TRUE" ]; then
    echo "Getting bin by barcode matrix as well ..."
    bin_file=${mat_dir}/${OUTPUT_PREFIX}_bin.bed
    mtx_bin_dir=${mat_dir}/bin_mat_resl${BIN_RESL}
    mkdir -p $mtx_bin_dir
    ${BEDTOOLS_PATH}/bedtools makewindows -g $CHROM_SIZE_FILE -w $BIN_RESL > $bin_file
    ${R_PATH}/R --vanilla --args ${OUTPUT_DIR}/summary/fragments.bed $bin_file ${mtx_bin_dir} < ${curr_dir}/src/get_mtx.R 
    rm $bin_file
fi


#rm ${mat_dir}/tmp.sam

echo "Get Matrix Done!"






