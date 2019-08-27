#!/bin/bash

set -e
unset PYTHONPATH

input_fastqs=$1

curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf $2
read_conf $3

## 2.trimming
echo "Trimming ..."
fastqs=(${input_fastqs//,/ })
dfastq1=$(basename ${fastqs[0]})
dfastq2=$(basename ${fastqs[1]})
${curr_dir}/trimming.sh ${fastqs[0]},${fastqs[1]} $2 $3


## 3.mapping 
echo "Mapping ..."
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
feature_file=${OUTPUT_DIR}/peaks/${PEAK_CALLER}/${OUTPUT_PREFIX}_features_BlacklistRemoved.bed
${curr_dir}/get_mtx.sh $feature_file $2 $3 &

echo "QC per cell ..."
frag_file=${OUTPUT_DIR}/summary/${OUTPUT_PREFIX}.fragments.bed
${curr_dir}/qc_per_barcode.sh $frag_file $2 $3 &

wait

## 7. call cell
echo "call cell ..."
mat_file=${OUTPUT_DIR}/raw_matrix/${PEAK_CALLER}/matrix.mtx
${curr_dir}/call_cell.sh $mat_file $2 $3

## report preprocessing QC
echo "generating report ..."
${curr_dir}/report.sh $OUTPUT_DIR/summary $2 $3

