#!/bin/bash


inputs=$1  ## seurat obj,and expeted doublet rate,
 
# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"


inputs=(${inputs//,/ })
seuratObj_atac_file=${inputs[0]}
seuratObj_rna_file=${inputs[1]}
gtf_file=${inputs[2]}

curr_dir=`dirname $0`

${R_PATH}/Rscript --vanilla ${curr_dir}/src/labelTransfer.R $seuratObj_atac_file $seuratObj_rna_file $GENOME_NAME $gtf_file

echo "Transfer Cell_Type label from scRNA-seq done!"

