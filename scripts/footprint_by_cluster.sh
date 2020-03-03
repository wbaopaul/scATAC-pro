#!/bin/bash

set -e

inputs=$1  ## two input bam files, the first must be something like *cluster_*.bam

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

if [[ ! -d "$HINT_PATH" ]]; then
    which rgt-hint > /dev/null
    if [ $? = "0" ]; then
        HINT_PATH=$(dirname `which rgt-hint`)
    else
        echo "HINT_PATH not provided or detected, please install rgt-hint!"
        exit
    fi
fi

down_dir=${OUTPUT_DIR}/downstream_analysis/${PEAK_CALLER}/${CELL_CALLER}
output_dir=${down_dir}/footprint
mkdir -p $output_dir

bams=(${inputs//,/ })
bam1=${bams[0]}
bam2=${bams[1]}

echo "construct regions by DA ..."
prefix1=${bam1/.bam/}
prefix1=$(basename $prefix1)
prefix2=${bam2/.bam/}
prefix2=$(basename $prefix2)
seurat_obj=${down_dir}/seurat_obj.rds

if [[ $prefix1 != *"cluster_"* ]]; then
    echo "the first bam file should be like cluster_*.bam"
    exit
fi
if [[ $prefix2 != *"cluster_"* ]]; then
    prefix2=rest
fi

prefix1=${prefix1/cluster_/}   ## absolute cluster name
prefix2=${prefix2/cluster_/}

input_peak=${down_dir}/differential_accessible_features_${prefix1}_vs_${prefix2}.txt

if [ ! -e "$input_peak" ]; then
    echo "do DA..."
    ${R_PATH}/Rscript --vanilla ${curr_dir}/src/runDA.R $seurat_obj ${down_dir} $prefix1 $prefix2 wilcox 
    echo "DA done succefully!"
fi

unset PYTHONPATH


output_dir=${output_dir}/${prefix1}_vs_${prefix2} ## make a subfolder to host the result
mkdir -p $output_dir



echo "predict footprint..."
${HINT_PATH}/rgt-hint footprinting --atac-seq --paired-end --organism=${GENOME_NAME} --output-location=$output_dir  \
        --output-prefix=${prefix1} $bam1 $input_peak &

${HINT_PATH}/rgt-hint footprinting --atac-seq --paired-end --organism=${GENOME_NAME} --output-location=$output_dir  \
        --output-prefix=${prefix2} $bam2 $input_peak &

wait

echo "overlap with motif annotation ... "
${HINT_PATH}/rgt-motifanalysis matching --organism=${GENOME_NAME} --input-files ${output_dir}/${prefix1}.bed ${output_dir}/${prefix2}.bed --output-location $output_dir

echo "Differential binding analysis ..."
${HINT_PATH}/rgt-hint differential --organism=${GENOME_NAME} --bc --nc 4 --mpbs-files=${output_dir}/${prefix1}_mpbs.bed,${output_dir}/${prefix2}_mpbs.bed --reads-files=$bam1,$bam2 --conditions=$prefix1,$prefix2 --output-location=${output_dir}

if [ -d "$TMP0" ]; then
    rm -r TMP0
fi

echo "Footprinting Analysis Done!"


