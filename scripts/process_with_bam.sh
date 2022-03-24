#!/bin/bash

set -e
unset PYTHONPATH

bam_file=$1

curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf $2
read_conf $3

mapRes_dir=${OUTPUT_DIR}/mapping_result
qc_dir=${OUTPUT_DIR}/summary

## 1.generate mapping qc
bash ${curr_dir}/bam2qc.sh $bam_file $2 $3

echo "Getting fragments (read pairs) per barcode file..."
${PERL_PATH}/perl ${curr_dir}/src/simply_bam2frags.pl --read_file ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.MAPQ${MAPQ}.bam \
        --output_file ${qc_dir}/${OUTPUT_PREFIX}.fragments.tsv --samtools_path $SAMTOOLS_PATH

${TABIX_PATH}/bgzip -f ${qc_dir}/${OUTPUT_PREFIX}.fragments.tsv.len

echo "sort and index fragment file ..."
# faster sort by R data.table
${R_PATH}/Rscript --vanilla ${curr_dir}/src/sort_frags.R ${qc_dir}/${OUTPUT_PREFIX}.fragments.tsv

#sort -k1,1 -k2,2n -T ${mapRes_dir}/tmp/  ${qc_dir}/${OUTPUT_PREFIX}.fragments.tsv > ${qc_dir}/${OUTPUT_PREFIX}.fragments.sorted.tsv
#mv ${qc_dir}/${OUTPUT_PREFIX}.fragments.sorted.tsv ${qc_dir}/${OUTPUT_PREFIX}.fragments.tsv
${TABIX_PATH}/bgzip -f ${qc_dir}/${OUTPUT_PREFIX}.fragments.tsv 
${TABIX_PATH}/tabix -f -p bed ${qc_dir}/${OUTPUT_PREFIX}.fragments.tsv.gz


## 2.call peak
echo "Calling peaks ..."
${curr_dir}/call_peak.sh $bam_file $2 $3 &

## 3.generate aggregated signal
echo "generating aggregated signal ..."
${curr_dir}/generate_signal.sh $bam_file $2 $3 &
wait

## 4.generate matrix
echo "generating raw matrix and qc per barcode..."
frag_file=${OUTPUT_DIR}/summary/${OUTPUT_PREFIX}.fragments.tsv.gz
feature_file=${OUTPUT_DIR}/peaks/${PEAK_CALLER}/${OUTPUT_PREFIX}_features_BlacklistRemoved.bed
${curr_dir}/get_mtx.sh ${frag_file},${feature_file} $2 $3 &

echo "QC per cell ..."
qc_inputs=${frag_file},${feature_file}
${curr_dir}/qc_per_barcode.sh $qc_inputs $2 $3 &

wait

## 5. call cell
echo "call cell ..."
mat_file=${OUTPUT_DIR}/raw_matrix/${PEAK_CALLER}/matrix.mtx
${curr_dir}/call_cell.sh $mat_file $2 $3

## 6. remove doublets
filtered_mtx_file=${OUTPUT_DIR}/filtered_matrix/${PEAK_CALLER}/${CELL_CALLER}/matrix.rds
${curr_dir}/rmDoublets.sh ${filtered_mtx_file},0.03 $2 $3

## 7. mapping qc for cell barcodes
map_dir=${OUTPUT_DIR}/mapping_result
input_bam=${map_dir}/${OUTPUT_PREFIX}.positionsort.bam
input_bc=${OUTPUT_DIR}/filtered_matrix/${PEAK_CALLER}/${CELL_CALLER}/barcodes_doubletsRemoved.txt
${curr_dir}/get_bam4Cells.sh ${input_bam},${input_bc} $2 $3

## 8.report preprocessing QC
echo "generating report ..."
${curr_dir}/report.sh ${OUTPUT_DIR}/summary $2 $3


