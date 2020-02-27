#!/bin/bash


mtx_files=$1  ## mtx for each sample, separated by ,
 
# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

output_dir=${OUTPUT_DIR}/integrated
mkdir -p $output_dir

curr_dir=`dirname $0`

${R_PATH}/Rscript --vanilla ${curr_dir}/src/integrate_seu.R $mtx_files $K_CLUSTERS $output_dir $GENOME_NAME $TSS $norm_by $REDUCTION $nREDUCTION $Top_Variable_Features $Integrate_By 

abs_out_dir=`cd ${output_dir}; pwd`

if [ "$prepCello4Integration" = "TRUE" ]; then
    seurat_file=${abs_out_dir}/seurat_obj_${Integrate_By}.rds
    assay4cello=ATAC
    if [[ $Integrate_By == "seurat"  ]]; then
        assay4cello=integrated
    fi
    ${R_PATH}/Rscript --vanilla ${curr_dir}/src/interface2cello.R $seurat_file $assay4cello
    ## write config file
    organism=hsa
    if [[ $GENOME_NAME =~ "mm" ]]; then
        organism=mmu
    fi

    echo "default:" > ${abs_out_dir}/VisCello_obj/config.yml
    echo "  study_name: $OUTPUT_PREFIX_integrated " >> ${abs_out_dir}/VisCello_obj/config.yml
    echo "  study_description: NNN " >> ${abs_out_dir}/VisCello_obj/config.yml
    echo "  organism: $organism " >> ${abs_out_dir}/VisCello_obj/config.yml
    echo "  feature_name_column: 'symbol' " >> ${abs_out_dir}/VisCello_obj/config.yml
    echo "  feature_id_column: 'symbol' " >> ${abs_out_dir}/VisCello_obj/config.yml
fi
