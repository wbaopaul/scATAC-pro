#!/bin/bash

## integrate from feature/peak files
set -e

input_peaks=$1  ## two or more input features files for two or more samples, seperated by ,

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

peaks=(${input_peaks//,/ })
peakLength=${#peaks[@]}

## put output into integrated_dir
integrated_dir=${OUTPUT_DIR}/integrated
mkdir -p ${integrated_dir}

peak_dir=${integrated_dir}/peaks
mkdir -p $peak_dir

reConst_mtx_dir=${integrated_dir}/reConstructed_matrix
mkdir -p $reConst_mtx_dir

echo "merge peaks ..."
feature_file=${peak_dir}/merged_peaks.bed
${R_PATH}/R --vanilla --args $input_peaks $feature_file < ${curr_dir}/src/mergePeaks.R

echo "ReConstructing peak-by-cell matrix for each sample ..."
echo "Using the previously called cells and merged peaks ..."
## supporse each sample was constructed by scATAC-pre
## so the fragment files are saved correspondly

ABS_PATH=`cd "$OUTPUT_DIR"; pwd`

## reconstruct the peak-cell matrix
## not re-call cells
mtx_files='TMP' 
re='^[0-9]+([.][0-9]+)?$'
for pk0 in "${peaks[@]}"
do 
    if [[ $pk0 =~ $re  ]] 2>/dev/null;then
        continue  ## the gap and qvalue parameter
    fi
    echo "Reconstruct matrix for sample related to $pk0: "
    sample0=$(basename $pk0)
    sample0=`echo $sample0 | awk -F. '{print $1}'`
    sample0=${sample0/_features_BlacklistRemoved/}
    pk0_dir=$(dirname $pk0)
    frag0_dir=`cd "$pk0_dir"; cd "../../summary"; pwd`       
    mat0_dir=`cd "$pk0_dir"; cd "../../filtered_matrix"; pwd`       
    #frag0_file=$(find $frag0_dir -name "*fragments.tsv.gz")
    frag0_file=$(find $frag0_dir -name "*fragments*" | grep -v "\.len" | grep -v tbi)
    mat0_dir=${mat0_dir}/${PEAK_CALLER}/${CELL_CALLER}
    #bc0_file=$(find ${mat0_dir} -name "*barcodes.txt")
    bc0_file=${mat0_dir}/barcodes_doubletsRemoved.txt
    if [ ! -e "$bc0_file" ]; then
         bc0_file=${mat0_dir}/barcodes.txt 
    fi
    
    new_mtx0_dir=${reConst_mtx_dir}/${sample0}
    mkdir -p $new_mtx0_dir
    bash ${curr_dir}/reConstMtx.sh ${feature_file},${frag0_file},${bc0_file},${new_mtx0_dir} $2 $3
    mtx_files=${mtx_files},${new_mtx0_dir}/matrix.rds
done

echo "Integrate by Seurat v3 ..."
echo -e "These are new mtx files: $mtx_files"
mtx_files=${mtx_files/TMP,/}

${R_PATH}/Rscript --vanilla ${curr_dir}/src/integrate_mtx.R $mtx_files $K_CLUSTERS $integrated_dir $GENOME_NAME $TSS $norm_by $REDUCTION $nREDUCTION $nFeature4Integration $Integrate_By

abs_out_dir=`cd ${integrated_dir}; pwd`

if [ "$prepCello4Integration" = "TRUE" ]; then
    seurat_file=${abs_out_dir}/seurat_obj_${Integrate_By}.rds
    assay4cello=ATAC
    if [[ $Integrate_By == "seurat"  ]]; then
        assay4cello=integrated
    fi
    ${R_PATH}/Rscript --vanilla ${curr_dir}/src/interface2cello.R $seurat_file $assay4cello $TSS
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

