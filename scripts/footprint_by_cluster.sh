#!/bin/bash

set -e

inputs=$1  ## two input bam file for two cluster, seperated by ,

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

output_dir=${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}/data_by_cluster/footprint
mkdir -p $output_dir

bams=(${inputs//,/ })

bam1=${bams[0]}
bam2=${bams[1]}

if [[ ! -d "$HINT_PATH" ]]; then
    which rgt-hint > /dev/null
    if [ $? = "0" ]; then
        HINT_PATH=$(dirname `which rgt-hint`)
    else
        echo "HINT_PATH not provided or detected, please install rgt-hint!"
        exit
    fi
fi


echo "construct regions by DA ..."
prefix1=${bam1/.bam/}
prefix1=$(basename $prefix1)
prefix2=${bam2/.bam/}
prefix2=$(basename $prefix2)
seurat_obj=${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}/seurat_obj.rds

cl1=${prefix1/cluster_/}   ## absolute cluster name
cl2=${prefix2/cluster_/}

if [[ $prefix1 == *"cluster_"* ]] && [[ $prefix2 == *"cluster_"* ]]; then
    input_peak=${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}/differential_accessible_features_cluster_${cl1}_VS_cluster_${cl2}.txt
    cl1_tmp=$cl1
    cl2_tmp=$cl2
else
    input_peak=${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}/differential_accessible_features_cluster_one_VS_cluster_rest.txt
    cl1_tmp=one
    cl2_tmp=rest
fi


if [ ! -e "$input_peak" ]; then
    echo "do DA..."
    ${R_PATH}/Rscript --vanilla ${curr_dir}/src/runDA.R $seurat_obj ${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER} $cl1_tmp $cl2_tmp wilcox 
    echo "DA done succefully!"
fi


unset PYTHONPATH

if [[ $prefix1 != *"cluster"* ]]; then
    prefix1=rest
fi

if [[ $prefix2 != *"cluster"* ]]; then
    prefix2=rest
fi

echo "predict footprint..."
${HINT_PATH}/rgt-hint footprinting --atac-seq --paired-end --organism=${GENOME_NAME} --output-location=$output_dir  \
        --output-prefix=${prefix1} $bam1 $input_peak &

${HINT_PATH}/rgt-hint footprinting --atac-seq --paired-end --organism=${GENOME_NAME} --output-location=$output_dir  \
        --output-prefix=${prefix2} $bam2 $input_peak &

wait

echo "overlap with motif annotation ... "
${HINT_PATH}/rgt-motifanalysis matching --organism=${GENOME_NAME} --input-files ${output_dir}/${prefix1}.bed ${output_dir}/${prefix2}.bed --output-location $output_dir

echo "Differential binding analysis ..."
${HINT_PATH}/rgt-hint differential --organism=${GENOME_NAME} --bc --nc 4 --mpbs-files=${output_dir}/${prefix1}_mpbs.bed,${output_dir}/${prefix2}_mpbs.bed --reads-files=$bam1,$bam2 --conditions=$prefix1,$prefix2 --output-location=${output_dir}/${prefix1}_${prefix2}

if [ -d "$TMP0" ]; then
    rm -r TMP0
fi

echo "Footprinting Analysis Done!"


