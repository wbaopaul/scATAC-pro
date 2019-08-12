#!/bin/bash

set -e

input_bam=$1  ## two ior more input bam files for two sample, seperated by ,

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

bams=(${input_bam//,/ })

bam1=${bams[0]}
bam2=${bams[1]}

output_dir=${OUTPUT_DIR}/downstream_analysis/${CELL_CALLER}/data_by_sample/footprint
mkdir -p $output_dir


echo "call peaks..."
unset PYTHONPATH
peaks_dir=${OUTPUT_DIR}/peaks/MACS2/integrated

for bam0 in "${bams[$@]}"; do
    mkdir -p $peaks_dir
    sample0=${bam0/.bam/}
    sample0=$(basename $sample0)
    ${MACS2_PATH}/macs2 callpeak -t $bam0 --outdir $peaks_dir -n $sample0 -f BAM $MACS2_OPTS

    ## remove peaks overlapped with blacklist
    ${BEDTOOLS_PATH}/bedtools intersect -a ${work_dir}/${sample0}_peaks.narrowPeak -b $BLACKLIST -v \
        > ${work_dir}/${sample0}_features_BlacklistRemoved.bed
done

echo "merge peaks ..."
${R_PATH}/R --vanilla --args $work_dir < ${curr_dir}/src/merge_peaks.R
mv ${work_dir}/merged_peaks.bed ${work_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}_features_BlacklistRemoved.bed


echo "Constructing raw peak-by-cell matrix for each sample ..."
## supporse bam was constructed by scATAC-pre
## so the fragment files are saved correspondly
feature_file=${work_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}_features_BlacklistRemoved.bed

ABS_PATH=`cd "$OUTPUT_DIR"; pwd`

for bam0 in "${bams[$@]}"; do 
    sample0=${bam0/.bam/}
    sample0=$(basename $sample0)
    raw_mtx_dir=${OUTPUT_DIR}/raw_matrix/${sample0}
    ${R_PATH}/R --vanilla --args ${OUTPUT_DIR}/summary/${sample0}.fragments.bed $feature_file ${raw_mtx_dir} 5 5 < ${curr_dir}/src/get_mtx.R

    echo "QC per barcode for each sample ..."=
    ${R_PATH}/R --vanilla --args ${OUTPUT_DIR}/summary/${sample0}.fragments.bed $feature_file $PROMOTERS $TSS $ENHANCERS ${qc_dir}/${sample0}.qc_per_barcode.txt < ${curr_dir}/src/get_qc_per_barcode.R

    echo "Call cell for each sample ..."
    bc_stat_file=${ABS_PATH}/summary/${sample0}.qc_per_barcode.bed
    mapping_qc_file=${ABS_PATH}/summary/${sample0}.MappingStats
    fragments_file=${ABS_PATH}/summary/${sample0}.fragments.bed
    output_dir=${ABS_PATH}/filtered_matrix/${sample0}

    echo "Integrate by Seurat v3 ..."
done




