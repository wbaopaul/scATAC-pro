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


curr_dir=`dirname $0`
work_dir=`pwd`



${R_PATH}/Rscript --vanilla ${curr_dir}/src/render2report.R \
    ${abs_report_dir}/scATAC-pro_report_${OUTPUT_PREFIX}.html  $abs_out_dir ${work_dir}/${2}
echo "Report generation Done!"
