#!/bin/bash

## integrate from feature/peak files, given SampleSheet.csv file
set -e

input_csv=$1 ## input cvs file 

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

## put output into integrated_dir
integrated_dir=${OUTPUT_DIR}/integrated
mkdir -p ${integrated_dir}

##merged peak directory
merged_peak_dir=${integrated_dir}/peaks
reConst_mtx_dir=${integrated_dir}/reConstructed_matrix
mkdir -p $merged_peak_dir
mkdir -p $reConst_mtx_dir

echo "merge peaks ..."
feature_file=${merged_peak_dir}/merged_peaks.bed
${R_PATH}/R --vanilla --args ${input_csv},${mergePeaksWithin},${filterPeaksQvalue} $feature_file < ${curr_dir}/src/mergePeaks_csv.R

echo "ReConstructing peak-by-cell matrix for each sample ..."
echo "Using given called cells and merged peaks ..."
## supporse each sample was constructed by scATAC-pro
## so the fragment files are saved correspondly

ABS_PATH=`cd "$OUTPUT_DIR"; pwd`

## prepare inputs for each sample
sample_names=($(cut -d , -f1 ${input_csv} | sed 1d))
peak_files=($(cut -d , -f2 ${input_csv}| sed 1d))
frag_files=($(cut -d , -f3 ${input_csv} | sed 1d))
barcodes_files=($(cut -d , -f4 ${input_csv} | sed 1d))
nSample=${#peak_files[@]}

## reconstruct the peak-cell matrix
## not re-call cells
mtx_files='TMP' 
for (( i=0; i<${nSample}; i++  ))
do 
    sample0=${sample_names[$i]}
    pk0_file=${peak_files[$i]}
    frag0_file=${frag_files[$i]}
    bc0_file=${barcodes_files[$i]}
    echo "Reconstruct matrix for sample related to ${pk0_file}: "
    mtx0_dir=${reConst_mtx_dir}/${sample0}
    mkdir -p $mtx0_dir

    ## if frag file not specified 
    if [[ -z "$frag0_file" ]]; then
        sample0=$(basename $pk0_file)
        sample0=`echo $sample0 | awk -F. '{print $1}'`
        sample0=${sample0/_features_BlacklistRemoved/}
        pk0_dir=$(dirname $pk0_file)
        frag0_dir=`cd "$pk0_dir"; cd "../../summary"; pwd`       
        frag0_file=$(find $frag0_dir -name "*fragments*" | grep -v "\.len" | grep -v tbi)
    fi
    
    ## if barcodes file not specified
    if [[ -z "$b0_file" ]]; then
        sample0=$(basename $pk0_file)
        sample0=`echo $sample0 | awk -F. '{print $1}'`
        sample0=${sample0/_features_BlacklistRemoved/}
        pk0_dir=$(dirname $pk0_file)
        bc0_dir=`cd "$pk0_dir"; cd "../../filtered_matrix"; pwd`       
        bc0_dir=${bc0_dir}/${PEAK_CALLER}/${CELL_CALLER}
        bc0_file=${bc0_dir}/barcodes_doubletsRemoved.txt
        if [ ! -e "$bc0_file" ]; then
             bc0_file=${bc0_dir}/barcodes.txt 
        fi
    fi

    bash ${curr_dir}/reConstMtx.sh ${feature_file},${frag0_file},${bc0_file},${mtx0_dir} $2 $3
    mtx_files=${mtx_files},${mtx0_dir}/matrix.rds
done

echo "Integrate by Seurat v3 ..."
echo -e "These are new mtx files: $mtx_files"
mtx_files=${mtx_files/TMP,/}

${R_PATH}/Rscript --vanilla ${curr_dir}/src/integrate_mtx_csv.R $mtx_files $K_CLUSTERS $integrated_dir $GENOME_NAME $TSS $norm_by $REDUCTION $nDim4Integration $nFeature4Integration $Integrate_By $input_csv

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

