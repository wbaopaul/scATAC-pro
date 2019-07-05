#!/bin/bash

set -e
unset PYTHONPATH

input_fastqs=$1
input_fastqs=${input_fastqs// /} ## remove space

curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf $2
read_conf $3

outfile_prefix=${OUTPUT_PREFIX}.${MAPPING_METHOD}

## 1.demultiplexing
make INPUT_FILE=$input_fastqs CONFIG_FILE=$2 CONFIG_SYS=$3 demplx_fastq > ${logDir}/demplx.log 2>&1  

## 2.trimming
fastqs=(${input_fastqs//,/ })
dfastq1=demplxed_$(basename ${fastqs[0]})
dfastq2=demplxed_$(basename ${fastqs[1]})
fq1=${OUTPUT_DIR}/demplxed_fastq/${dfastq1}
fq2=${OUTPUT_DIR}/demplxed_fastq/${dfastq2}
make INPUT_FILE=${fq1},${fq2} CONFIG_FILE=$2 CONFIG_SYS=$3 trimming > ${logDir}/trimming.log 2>&1  

## 3.mapping 
if [ "$TRIM_METHOD" = "trim_galore" ]; then
    dfastq1_pre=`echo $dfastq1 | awk -F. '{print $1}'`
    dfastq2_pre=`echo $dfastq2 | awk -F. '{print $1}'`
    trimmed_fq1=${OUTPUT_DIR}/trimmed_fastq/${dfastq1_pre}_val_1.fq.gz
    trimmed_fq2=${OUTPUT_DIR}/trimmed_fastq/${dfastq1_pre}_val_2.fq.gz
    mapping_inputs=${trimmed_fq1},${trimmed_fq2}
elif [ "$TRIM_METHOD" = "Trimmomatic" ]; then
    trimmed_fq1=${OUTPUT_DIR}/trimmed_fastq/trimmed_paired_${dfastq1}
    trimmed_fq2=${OUTPUT_DIR}/trimmed_fastq/trimmed_paired_${dfastq2}
    mapping_inputs=${trimmed_fq1},${trimmed_fq2}
else
    mapping_inputs=${fq1},${fq2}
fi
 
make INPUT_FILE=$mapping_inputs CONFIG_FILE=$2 CONFIG_SYS=$3 mapping > ${logDir}/mapping.log 2>&1  

## 4.call peak
bam_file=${OUTPUT_DIR}/mapping_result/${outfile_prefix}.positionsort.MAPQ${MAPQ}.bam
make INPUT_FILE=$bam_file CONFIG_FILE=$2 CONFIG_SYS=$3 call_peak > ${logDir}/call_peak.log 2>&1  

## 5.generate matrix
if [ ! $USE_BIN ]; then
    feature_file=${OUTPUT_DIR}/peaks/${PEAK_CALLER}/${outfile_prefix}_peaks_BlacklistRemoved.bed
else
    feature_file="xx"
fi
make INPUT_FILE=$feature_file CONFIG_FILE=$2 CONFIG_SYS=$3 get_mtx > ${logDir}/get_mtx.log 2>&1 & 

## 6.generate aggregated signal
make INPUT_FILE=$bam_file CONFIG_FILE=$2 CONFIG_SYS=$3 aggr_signal > ${logDir}/aggr_signal.log 2>&1 & 

wait


## qc and call cell
if [ $USE_BIN ]; then
    mat_file=${OUTPUT_DIR}/raw_matrix/bin_mat_resl${BIN_RESL}/${outfile_prefix}.bin_resl${BIN_RESL}.mtx
else
    mtt_file=${OUTPUT_DIR}/raw_matrix/peak_mat/${outfile_prefix}.peak.barcode.mtx
fi

frag_file=${OUTPUT_DIR}/raw_matrix/fragments.bed

if [ "$CELL_CALLER" != "filtering" ]; then
    ## 7.qc per barcode
    make INPUT_FILE=$frag_file CONFIG_FILE=$2 CONFIG_SYS=$3 qc_per_barcode > ${logDir}/qc_per_barcode.log 2>&1 &

    ## 8.call cell
    make INPUT_FILE=$mat_file CONFIG_FILE=$2 CONFIG_SYS=$3 call_cell > ${logDir}/call_cell.log 2>&1 &
else
    ## 7.qc per barcode
    make INPUT_FILE=$frag_file CONFIG_FILE=$2 CONFIG_SYS=$3 qc_per_barcode > ${logDir}/qc_per_barcode.log 2>&1 

    ## 8.call cell
    make INPUT_FILE=$mat_file CONFIG_FILE=$2 CONFIG_SYS=$3 call_cell > ${logDir}/call_cell.log 2>&1

fi

## report preprocessing QC


