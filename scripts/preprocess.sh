#!/bin/bash

set -e
unset PYTHONPATH

input_fastqs=$1

curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf $2
read_conf $3

outfile_prefix=${OUTPUT_PREFIX}.${MAPPING_METHOD}

## 1.demultiplexing
#make INPUT_FILE=$input_fastqs CONFIG_FILE=$2 CONFIG_SYS=$3 demplx_fastq > ${logDir}/demplx.log 2>&1  
${curr_dir}/dex_fastq.sh $1 $2 $3

## 2.trimming
fastqs=(${input_fastqs//,/ })
dfastq1=demplxed_$(basename ${fastqs[0]})
dfastq2=demplxed_$(basename ${fastqs[1]})
fq1=${OUTPUT_DIR}/demplxed_fastq/${dfastq1}
fq2=${OUTPUT_DIR}/demplxed_fastq/${dfastq2}
#make INPUT_FILE=${fq1},${fq2} CONFIG_FILE=$2 CONFIG_SYS=$3 trimming > ${logDir}/trimming.log 2>&1  
${curr_dir}/trimming.sh ${fq1},${fq2} $2 $3

## 3.mapping 
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
    mapping_inputs=${fq1},${fq2}
fi
 
${curr_dir}/mapping.sh $mapping_inputs $2 $3

## 4.call peak
bam_file=${OUTPUT_DIR}/mapping_result/${outfile_prefix}.positionsort.MAPQ${MAPQ}.noDuplicates.bam
${curr_dir}/call_peak.sh $bam_file $2 $3

## 5.generate matrix

feature_file=${OUTPUT_DIR}/peaks/${PEAK_CALLER}/${outfile_prefix}_features_BlacklistRemoved.bed

${curr_dir}/get_mtx.sh $feature_file $2 $3 &

## 6.generate aggregated signal
${curr_dir}/generate_signal.sh $bam_file $2 $3 &
wait


## qc and call cell
mat_file=${OUTPUT_DIR}/raw_matrix/${PEAK_CALLER}/matrix.mtx

frag_file=${OUTPUT_DIR}/summary/fragments.bed

if [ "$CELL_CALLER" != "filtering" ]; then
    ## 7.qc per barcode
    ${curr_dir}/qc_per_barcode.sh $frag_file $2 $3 &
    ## 8.call cell

    ${curr_dir}/call_cell.sh $mat_file $2 $3 &
else
    ## 7.qc per barcode
    
    ${curr_dir}/qc_per_barcode.sh $frag_file $2 $3 
    ## 8.call cell

    ${curr_dir}/call_cell.sh $mat_file $2 $3 
fi

wait
## report preprocessing QC
${curr_dir}/report.sh $OUTPUT_DIR/summary $2 $3
