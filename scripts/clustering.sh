#!/bin/bash


mtx_file=$1
 
# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

output_dir=${OUTPUT_DIR}/downstream_analysis/${CELL_CALLER}
mkdir -p $output_dir

curr_dir=`dirname $0`

${R_PATH}/Rscript --vanilla ${curr_dir}/src/clustering.R $mtx_file $CLUSTERING_METHOD $K_CLUSTERS $output_dir $GENOME_NAME $TSS $norm_by $REDUCTION $nREDUCTION $Top_Variable_Features 

if [ "$prepCello" = "TRUE" ]; then
    seurat_file=${output_dir}/seurat_obj.rds
    ${R_PATH}/Rscript --vanilla ${curr_dir}/src/interface2cello.R $seurat_file
    ## write config file
    organism=hsa
    if [[ $GENOME_NAME =~ "mm" ]]; then
        organism=mmu
    fi

    echo "default:" > ${output_dir}/VisCello_obj/config.yml
    echo "  study_name: $OUTPUT_PREFIX " >> ${output_dir}/VisCello_obj/config.yml
    echo "  study_description: NNN " >> ${output_dir}/VisCello_obj/config.yml
    echo "  organism: $organism " >> ${output_dir}/VisCello_obj/config.yml
    echo "  feature_name_column: 'symbol' " >> ${output_dir}/VisCello_obj/config.yml
    echo "  feature_id_column: 'symbol' " >> ${output_dir}/VisCello_obj/config.yml
fi
