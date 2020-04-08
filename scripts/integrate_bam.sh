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

#OUTPUT_DIR=${OUTPUT_DIR}/integrated
#mkdir -p $OUTPUT_DIR
echo "call peaks..."
unset PYTHONHOME
unset PYTHONPATH
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
${R_PATH}/R --vanilla --args $peak_dir < ${curr_dir}/src/merge_peaks_cls.R


## remove peaks overlapped with blacklist
${BEDTOOLS_PATH}/bedtools intersect -a ${peak_dir}/merged_peaks.bed -b $BLACKLIST -v \
    > ${peak_dir}/${OUTPUT_PREFIX}_features_BlacklistRemoved.bed 

find ${peak_dir}/* -not -name "*Blacklist*" | grep -v narrow | xargs rm


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
    frag0_file=$(find $frag0_dir -name "*fragments.txt")
    ${R_PATH}/R --vanilla --args $frag0_file $feature_file ${raw_mtx_dir} 5 5 < ${curr_dir}/src/get_mtx.R &

    echo "QC per barcode for each sample ..."
    ${R_PATH}/R --vanilla --args $frag0_file $feature_file $PROMOTERS $TSS $ENHANCERS ${qc_dir}/${sample0}.MACS2.qc_per_barcode.txt < ${curr_dir}/src/get_qc_per_barcode.R &
    wait
done


# call cell
mtx_files="TMP"
for bam0 in "${bams[@]}"
do 
    sample0=$(basename $bam0)
    sample0=`echo $sample0 | awk -F. '{print $1}'`
    bam0_dir=$(dirname $bam0)
    frag0_dir=`cd "$bam0_dir"; cd "../summary"; pwd`       
    echo "Call cell for each sample ..."
    bc_stat_file=${ABS_PATH}/summary/${sample0}.MACS2.qc_per_barcode.txt
    filtered_mtx_dir=${ABS_PATH}/filtered_matrix/${PEAK_CALLER}/${CELL_CALLER}/${sample0}
    raw_mtx_dir=${ABS_PATH}/raw_matrix/${sample0}
    input_mtx=${raw_mtx_dir}/matrix.mtx

    mkdir -p $filtered_mtx_dir
    if [ ${CELL_CALLER} = 'EmptyDrop' ];then
        ${R_PATH}/R --vanilla --args $input_mtx $filtered_mtx_dir ${EmptyDrop_FDR} < ${curr_dir}/src/EmptyDrop.R &
    fi
    
    if [ ${CELL_CALLER} = 'FILTER' ];then

        echo "${R_PATH}/Rscript  ${curr_dir}/src/filter_barcodes.R  --bc_stat_file $bc_stat_file --raw_mtx_file $input_mtx --output_dir $filtered_mtx_dir ${FILTER_BC_CUTOFF}" > tmpJob_${sample0}
        bash tmpJob_${sample0} &
    fi
    mtx_files=${mtx_files},${filtered_mtx_dir}/matrix.mtx
done
wait
rm tmpJob*

echo "Integrate by Seurat v3 ..."
mtx_files=${mtx_files/TMP,/}
#${curr_dir}/integrate_seu.sh $mtx_files $2 $3

${R_PATH}/Rscript --vanilla ${curr_dir}/src/integrate_seu.R $mtx_files $CLUSTERING_METHOD $K_CLUSTERS $OUTPUT_DIR $GENOME_NAME $TSS

