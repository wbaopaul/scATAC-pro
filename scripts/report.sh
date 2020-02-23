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

#grep = $2 | grep -v ^# | awk -F= '{print $1}' | awk '{$1=$1;print}' > ${abs_out_dir}/vrs.txt 

#grep = $2 | grep -v ^# | awk -F= '{print $2}' | awk -F# '{print $1}' | awk '{$1=$1;print}' > ${abs_out_dir}/vls.txt


${R_PATH}/Rscript --vanilla ${curr_dir}/src/render2report.R \
    ${abs_report_dir}/scATAC-pro_report_${OUTPUT_PREFIX}.html  $abs_out_dir ${work_dir}/${2}
rm ${abs_out_dir}/vrs.txt
rm ${abs_out_dir}/vls.txt         
echo "Report generation Done!"
