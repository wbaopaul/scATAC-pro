#!/bin/bash


mtx_file=$1
 
# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

abs_out_dir=`cd ${OUTPUT_DIR}; pwd`
abs_down_dir=${abs_out_dir}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}
mkdir -p $abs_down_dir

curr_dir=`dirname $0`

${R_PATH}/Rscript --vanilla ${curr_dir}/src/clustering.R $mtx_file $CLUSTERING_METHOD $K_CLUSTERS $abs_down_dir $GENOME_NAME $TSS $norm_by $REDUCTION $nREDUCTION $Top_Variable_Features 


if [ "$prepCello" = "TRUE" ]; then
    seurat_file=${abs_down_dir}/seurat_obj.rds
    ${R_PATH}/Rscript --vanilla ${curr_dir}/src/interface2cello.R $seurat_file ATAC
    ## write config file
    organism=hsa
    if [[ $GENOME_NAME =~ "mm" ]]; then
        organism=mmu
    fi

    echo "default:" > ${abs_down_dir}/VisCello_obj/config.yml
    echo "  study_name: $OUTPUT_PREFIX " >> ${abs_down_dir}/VisCello_obj/config.yml
    echo "  study_description: NNN " >> ${abs_down_dir}/VisCello_obj/config.yml
    echo "  organism: $organism " >> ${abs_down_dir}/VisCello_obj/config.yml
    echo "  feature_name_column: 'symbol' " >> ${abs_down_dir}/VisCello_obj/config.yml
    echo "  feature_id_column: 'symbol' " >> ${abs_down_dir}/VisCello_obj/config.yml
   
fi
