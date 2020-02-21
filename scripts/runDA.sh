#!/bin/bash


groups=$1  ## group1,group2, like 0;1,2 which compare cluster 0,1 with cluster2
 
# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

output_dir=${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}
seuratObj_file=${output_dir}/seurat_obj.rds
curr_dir=`dirname $0`

groups=(${groups//,/ })
group1=${groups[0]}
group2=${groups[1]}
${R_PATH}/Rscript --vanilla ${curr_dir}/src/runDA.R $seuratObj_file $output_dir $group1 $group2 $test_use

echo "Differential analysis done!"
