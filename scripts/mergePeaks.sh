#!/bin/bash

set -e

## merge multiple peaks
##  
## given all peak files
## and a distance parameter like 200, separated by comma
## the merged peaks will be saved in output/peaks/merged_peaks.bed
## example inputs: peak_File1,peak_File2,peak_File3,200
## output file: output/peaks/merged_peaks.bed

inputs=$1  

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"
mkdir -p ${OUTPUT_DIR}/peaks
merged_file="${OUTPUT_DIR}/peaks/merged_peaks.bed"

echo "merge peaks ..."
${R_PATH}/R --vanilla --args $inputs $merged_file < ${curr_dir}/src/mergePeaks.R

