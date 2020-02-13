#!/bin/bash


de_file=$1 ## the file provides forgroud genes for each cluster, all genes as background genes
## at least column 'gene_name' and 'cluster' should be provided
 
# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

output_dir=${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}
mkdir -p $output_dir

curr_dir=`dirname $0`

${R_PATH}/Rscript --vanilla ${curr_dir}/src/runGO.R $de_file $output_dir $GENOME_NAME $GO_TYPE

echo "GO analysis done!"
