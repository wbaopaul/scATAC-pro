#!/bin/bash

set -e
unset PYTHONPATH

input_fastqs=$1

curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf $2
read_conf $3

fastqs=(${input_fastqs//,/ })
## 2.trimming
echo "Trimming ..."
isSingleEnd=${isSingleEnd^^}
if [[ "$isSingleEnd"="FALSE"  ]]; then
    dfastq1=$(basename ${fastqs[0]})
    dfastq2=$(basename ${fastqs[1]})
    ${curr_dir}/trimming.sh ${fastqs[0]},${fastqs[1]} $2 $3
    if [ "$TRIM_METHOD" = "trim_galore" ]; then
        dfastq1_pre=`echo $dfastq1 | awk -F. '{print $1}'`
        dfastq2_pre=`echo $dfastq2 | awk -F. '{print $1}'`
        trimmed_fq1=${OUTPUT_DIR}/trimmed_fastq/${dfastq1_pre}_val_1.fq.gz
        trimmed_fq2=${OUTPUT_DIR}/trimmed_fastq/${dfastq2_pre}_val_2.fq.gz
        mapping_inputs=${trimmed_fq1},${trimmed_fq2}
    elif [ "$TRIM_METHOD" = "Trimmomatic" ]; then
        trimmed_fq1=${OUTPUT_DIR}/trimmed_fastq/trimmed_paired_${dfastq1}
        trimmed_fq2=${OUTPUT_DIR}/trimmed_fastq/trimmed_paired_${dfastq2}
        mapping_inputs=${trimmed_fq1},${trimmed_fq2}
    else
        mapping_inputs=${fastqs[0]},${fastqs[1]}
    fi
else
    dfastq1=$(basename ${fastqs[0]})
    ${curr_dir}/trimming.sh ${fastqs[0]} $2 $3
    if [ "$TRIM_METHOD" = "trim_galore" ]; then
        dfastq1_pre=`echo $dfastq1 | awk -F. '{print $1}'`
        trimmed_fq1=${OUTPUT_DIR}/trimmed_fastq/${dfastq1_pre}_trimmed.fq.gz
        mapping_inputs=${trimmed_fq1}
    elif [ "$TRIM_METHOD" = "Trimmomatic" ]; then
        trimmed_fq1=${OUTPUT_DIR}/trimmed_fastq/trimmed_${dfastq1}
        mapping_inputs=${trimmed_fq1}
    else
        mapping_inputs=${fastqs[0]}
    fi
fi

## 3.mapping 
echo "Mapping ..."

echo "Start Mapping ..." 
${curr_dir}/mapping.sh $mapping_inputs $2 $3

## 4.call peak
echo "Calling peaks ..."
bam_file=${OUTPUT_DIR}/mapping_result/${OUTPUT_PREFIX}.positionsort.MAPQ${MAPQ}.bam
${curr_dir}/call_peak.sh $bam_file $2 $3 &

## 5.generate aggregated signal
echo "generating aggregated signal ..."
${curr_dir}/generate_signal.sh $bam_file $2 $3 &
wait

## 6.generate matrix
echo "generating raw matrix and qc per barcode..."
frag_file=${OUTPUT_DIR}/summary/${OUTPUT_PREFIX}.fragments.txt
feature_file=${OUTPUT_DIR}/peaks/${PEAK_CALLER}/${OUTPUT_PREFIX}_features_BlacklistRemoved.bed
${curr_dir}/get_mtx.sh ${frag_file},${feature_file} $2 $3 &

echo "QC per cell ..."
qc_inputs=${frag_file},${feature_file}
${curr_dir}/qc_per_barcode.sh $qc_inputs $2 $3 &

wait

## 7. call cell
echo "call cell ..."
mat_file=${OUTPUT_DIR}/raw_matrix/${PEAK_CALLER}/matrix.mtx
${curr_dir}/call_cell.sh $mat_file $2 $3

## 8. mapping qc for cell barcodes
map_dir=${OUTPUT_DIR}/mapping_result
input_bam=${map_dir}/${OUTPUT_PREFIX}.positionsort.bam
input_bc=${OUTPUT_DIR}/filtered_matrix/${PEAK_CALLER}/${CELL_CALLER}/barcodes.txt
${curr_dir}/get_bam4Cells.sh ${input_bam},${input_bc} $2 $3

## report preprocessing QC
echo "generating report ..."
${curr_dir}/report.sh $OUTPUT_DIR/summary $2 $3

