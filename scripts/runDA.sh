#!/bin/bash


inputs=$1  ## seurat obj,group1,group2, like seurat_obj.rds,0:1,2 which compare cluster 0,1 with cluster2 in the seurat obj
 
# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

#output_dir=${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}
#oseuratObj_file=${output_dir}/seurat_obj.rds

inputs=(${inputs//,/ })
seuratObj_file=${inputs[0]}
group1=${inputs[1]}
group2=${inputs[2]}

output_dir=`dirname $seuratObj_file`
curr_dir=`dirname $0`

${R_PATH}/Rscript --vanilla ${curr_dir}/src/runDA.R $seuratObj_file $output_dir $group1 $group2 $test_use $TSS

echo "Differential analysis done!"
