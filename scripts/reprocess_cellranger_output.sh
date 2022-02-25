#!/bin/bash

set -e
unset PYTHONPATH

input_file=$1
files=(${input_file//,/ })
fileLength=${#files[@]}

if [ "$fileLength" != 2 ];then
    echo "Please specify a cellranger-atac generated bam file and a fragment file, separated by a comma!" >&2
    exit
fi

bam_file0=${files[0]}
frag_file0=${files[1]}

curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf $2
read_conf $3

mapRes_dir=${OUTPUT_DIR}/mapping_result
qc_dir=${OUTPUT_DIR}/summary

## 1.generate mapping qc
bash ${curr_dir}/bam2qc.sh $bam_file0 $2 $3

## 2.call peak
echo "Calling peaks ..."
bam_file=${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.bam
bam_file30=${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.MAPQ30.bam
${curr_dir}/call_peak.sh $bam_file30 $2 $3 &

## 3.generate aggregated signal
echo "generating aggregated signal ..."
${curr_dir}/generate_signal.sh $bam_file30 $2 $3 &

wait

## 4.generate matrix
echo "generating raw matrix and qc per barcode..."
frag_file=${OUTPUT_DIR}/summary/${OUTPUT_PREFIX}.fragments.tsv.gz
feature_file=${OUTPUT_DIR}/peaks/${PEAK_CALLER}/${OUTPUT_PREFIX}_features_BlacklistRemoved.bed
if [[ ! -f $frag_file ]]; then
    ln -s $frag_file0 $frag_file
fi
## softlink index file
if [[ -f ${frag_file0}.tbi ]]; then
    ln -s ${frag_file0}.tbi ${frag_file}.tbi
fi
${curr_dir}/get_mtx.sh ${frag_file},${feature_file} $2 $3 &

echo "QC per cell ..."
qc_inputs=${frag_file},${feature_file}
${curr_dir}/qc_per_barcode.sh $qc_inputs $2 $3 &

if [[ $CELL_MAP_QC = TRUE ]]; then
    echo "convert 10x bam into scATAC-pro bam ..."
    ${curr_dir}/convert10xbam_short.sh $bam_file,${bam_file}.tmp &
fi
wait

if [[ -f ${bam_file}.tmp ]]; then
    mv ${bam_file}.tmp $bam_file
    mv ${bam_file}.tmp.bai ${bam_file}.bai
fi

## 5. call cell
echo "call cell ..."
mtx_file=${OUTPUT_DIR}/raw_matrix/${PEAK_CALLER}/matrix.mtx
${curr_dir}/call_cell.sh $mtx_file $2 $3

## 6. remove doublets
input_bc=${OUTPUT_DIR}/filtered_matrix/${PEAK_CALLER}/${CELL_CALLER}/barcodes.txt
mtx_file=${OUTPUT_DIR}/filtered_matrix/${PEAK_CALLER}/${CELL_CALLER}/matrix.rds
if [[ $rmDoublets = TRUE  ]]; then
    ${curr_dir}/rmDoublets.sh ${mtx_file},0.03 $2 $3 
    input_bc=${OUTPUT_DIR}/filtered_matrix/${PEAK_CALLER}/${CELL_CALLER}/barcodes_doubletsRemoved.txt
    mtx_file=${OUTPUT_DIR}/filtered_matrix/${PEAK_CALLER}/${CELL_CALLER}/matrix_doubletsRemoved.rds
fi

## 7. mapping qc for cell barcodes 
if [[ $CELL_MAP_QC = TRUE ]]; then
   ## run get_bam4Cells
    ${curr_dir}/get_bam4Cells.sh ${bam_file},${input_bc} $2 $3 &
fi

## 8.clustering
${curr_dir}/clustering.sh ${mtx_file} $2 $3 & 

wait

## 9.report preprocessing QC
echo "generating report ..."
${curr_dir}/report.sh ${OUTPUT_DIR}/summary $2 $3

