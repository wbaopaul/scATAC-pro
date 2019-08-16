#!/bin/bash

set -e

input_bam=$1  ## two or more input bam files for two sample, seperated by ,

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

bams=(${input_bam//,/ })

bam1=${bams[0]}
bam2=${bams[1]}

OUTPUT_DIR=${OUTPUT_DIR}/integrated
mkdir -p $OUTPUT_DIR
echo "call peaks..."
unset PYTHONPAT
peak_dir=${OUTPUT_DIR}/peaks
mkdir -p $peak_dir


for bam0 in "${bams[@]}"
do
    sample0=$(basename $bam0)
    sample0=`echo $sample0 | awk -F. '{print $1}'`
    ${MACS2_PATH}/macs2 callpeak -t $bam0 --outdir $peak_dir -n $sample0 -f BAM $MACS2_OPTS & 
done

wait


echo "merge peaks ..."
${R_PATH}/R --vanilla --args $peak_dir < ${curr_dir}/src/merge_peaks.R


## remove peaks overlapped with blacklist
${BEDTOOLS_PATH}/bedtools intersect -a ${peak_dir}/merged_peaks.bed -b $BLACKLIST -v \
    > ${peak_dir}/${OUTPUT_PREFIX}_features_BlacklistRemoved.bed 
find ${peak_dir}/* -not -name "*Blacklist*" | xargs rm


echo "Constructing raw peak-by-cell matrix for each sample ..."
## supporse bam was constructed by scATAC-pre
## so the fragment files are saved correspondly
feature_file=${peak_dir}/${OUTPUT_PREFIX}_features_BlacklistRemoved.bed

ABS_PATH=`cd "$OUTPUT_DIR"; pwd`

## get raw matrix and qc
for bam0 in "${bams[@]}"
do 
    sample0=$(basename $bam0)
    sample0=`echo $sample0 | awk -F. '{print $1}'`
    raw_mtx_dir=${OUTPUT_DIR}/raw_matrix/${sample0}
    mkdir -p $raw_mtx_dir
    bam0_dir=$(dirname $bam0)
    frag0_dir=`cd "$bam0_dir"; cd "../summary"; pwd`       
    qc_dir=${ABS_PATH}/summary 
    mkdir -p $qc_dir
    frag0_file=$(find $frag0_dir -name "*fragments.bed")
    ${R_PATH}/R --vanilla --args $frag0_file $feature_file ${raw_mtx_dir} 5 5 < ${curr_dir}/src/get_mtx.R &

    echo "QC per barcode for each sample ..."
    ${R_PATH}/R --vanilla --args $frag0_file $feature_file $PROMOTERS $TSS $ENHANCERS ${qc_dir}/${sample0}.qc_per_barcode.txt < ${curr_dir}/src/get_qc_per_barcode.R &
    wait
done


# call cell
for bam0 in "${bams[@]}"
do 
    sample0=$(basename $bam0)
    sample0=`echo $sample0 | awk -F. '{print $1}'`
    bam0_dir=$(dirname $bam0)
    frag0_dir=`cd "$bam0_dir"; cd "../summary"; pwd`       
    echo "Call cell for each sample ..."
    bc_stat_file=${ABS_PATH}/summary/${sample0}.qc_per_barcode.txt
    output_dir=${ABS_PATH}/filtered_matrix/${sample0}
    input_mtx=${raw_mtx_dir}/matrix.mtx

    if [ ${CELL_CALLER} = 'EmptyDrop' ];then
        output_dir1=${output_dir}/EmptyDrop
        mkdir -p $output_dir1
        ${R_PATH}/R --vanilla --args $input_mtx $output_dir1 ${EmptyDrop_FDR} < ${curr_dir}/src/EmptyDrop.R &
    fi
    
    if [ ${CELL_CALLER} = 'FILTER' ];then
        output_dir1=${output_dir}/FILTER
        mkdir -p $output_dir1

        echo "${R_PATH}/Rscript  ${curr_dir}/src/filter_barcodes.R  --bc_stat_file $bc_stat_file --raw_mtx_file $input_mtx --output_dir $output_dir1 ${FILTER_BC_CUTOFF}" > tmpJob_${sample0}
        bash tmpJob_${sample0} &
    fi
done
wait
rm tmpJob*

echo "Integrate by Seurat v3 ..."

