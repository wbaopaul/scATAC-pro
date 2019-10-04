#!/bin/bash

set -e

input_bam=$1  ## two input bam file for two cluster, seperated by ,

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

bams=(${input_bam//,/ })

bam1=${bams[0]}
bam2=${bams[1]}

output_dir=${OUTPUT_DIR}/downstream_analysis/${CELL_CALLER}/data_by_cluster/footprint
mkdir -p $output_dir


if [[ ! -d "$HINT_PATH" ]]; then
    which rgt-hint > /dev/null
    if [ $? = "0" ]; then
        HINT_PATH=$(dirname `which rgt-hint`)
    else
        echo "HINT_PATH not provided or detected, please install rgt-hint!"
        exit
    fi
fi


## select regions (peaks) by DA ##
cluster1=${bam1/.bam/}
cluster1=$(basename $cluster1)
cluster2=${bam2/.bam/}
cluster2=$(basename $cluster2)
seurat_obj=${OUTPUT_DIR}/downstream_analysis/${CELL_CALLER}/seurat_obj.rds

cl1=${cluster1/cluster_/}
cl2=${cluster2/cluster_/}


input_peak=${OUTPUT_DIR}/downstream_analysis/${CELL_CALLER}/differential_peak_cluster_table.txt
if [ ! -e "$input_peak" ]; then
    echo "do DA..."
    mkdir -p TMP0
    ${R_PATH}/Rscript --vanilla ${curr_dir}/src/runDA.R $seurat_obj TMP0 all others LR 
    echo "DA done succefully!"
#    echo 'library(data.table); dd = fread("TMP0/differential_peak_cluster_table.txt");
#                dd[, "chr" := unlist(strsplit(peak, "-"))[1], by = peak];
#                dd[, "start" := as.integer(unlist(strsplit(peak, "-"))[2]), by = peak];
#                dd[, "end" := as.integer(unlist(strsplit(peak, "-"))[3]), by = peak];
#                dd = subset(dd, select = c("chr", "start", "end"));
#                setkey(dd, chr, start)
#                write.table(dd, file = "TMP0/select_peaks.bed", row.names = F, col.names = F,
#                            quote = F, sep = "\t") 
#    ' > TMP0/tmp.R
#    ${R_PATH}/R --vanilla < TMP0/tmp.R
    input_peak=TMP0/differential_peak_cluster_table.txt
fi


unset PYTHONPATH

echo "predict footprint..."
${HINT_PATH}/rgt-hint footprinting --atac-seq --paired-end --organism=${GENOME_NAME} --output-location=$output_dir  \
        --output-prefix=${cluster1} $bam1 $input_peak &

${HINT_PATH}/rgt-hint footprinting --atac-seq --paired-end --organism=${GENOME_NAME} --output-location=$output_dir  \
        --output-prefix=${cluster2} $bam2 $input_peak &

wait

echo "overlap with motif annotation ... "
${HINT_PATH}/rgt-motifanalysis matching --organism=${GENOME_NAME} --input-files ${output_dir}/${cluster1}.bed ${output_dir}/${cluster2}.bed --output-location $output_dir

echo "Differential binding analysis ..."
${HINT_PATH}/rgt-hint differential --organism=${GENOME_NAME} --bc --nc 4 --mpbs-file1=${output_dir}/${cluster1}_mpbs.bed \
    --mpbs-file2=${output_dir}/${cluster1}_mpbs.bed --reads-file1=$bam1 --reads-file2=$bam2 --condition1=$cluster1 --condition2=$cluster2 --output-location=${output_dir}/${cluster1}_${cluster2}

if [ -d "$TMP0" ]; then
    rm -r TMP0
fi

echo "Footprinting Analysis Done!"


