#!/bin/bash

### will search peak file in peaks/*narrowPeak under $2
set -e

input_bam=$1
region_file=${2}/peaks/${OUTPUT_PREFIX}.${MAPPING_METHOD}_peaks_filterBlacklist.narrowPeak
mat_dir="${2}/raw_matrix"
mkdir -p $mat_dir

ncore=$(nproc --all)
ncore=$(($ncore - 1))
## conver bam to sam
if [[ -z ${mat_dir}/tmp.sam ]]; then
  ${SAMTOOLS_PATH}/samtools view -h -@ $ncore  $input_bam > ${mat_dir}/tmp.sam
fi

echo "filter barcodes with less than 10 reads"
curr_dir=`dirname $0`
if [[ -z ${mat_dir}/tmp.bed ]]; then
	perl ${curr_dir}/cal_frac_mito.pl --read_file ${mat_dir}/tmp.sam  --read_length 50  --output_file ${mat_dir}/tmp.bed
	sed -i '1d' ${mat_dir}/tmp.bed
fi

awk '$2 > 10 {print $1}' ${mat_dir}/tmp.bed > ${mat_dir}/barcodes.reads.GT10

perl ${curr_dir}/extract_sam_gbarcodes.pl --barcode_file ${mat_dir}/barcodes.reads.GT10 --read_file ${mat_dir}/tmp.sam --output_file ${mat_dir}/barcodes.reads.GT10.sam

## use perl script to get the matrix
echo "getting matrix..."
${PERL_PATH}/perl ${curr_dir}/get_peak_barcode_mat.pl --region_file $region_file --read_file ${mat_dir}/barcodes.reads.GT10.sam --read_length 50 --output_file ${mat_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}.peak.barcode.mat  

rm ${mat_dir}/tmp.sam

echo "Get Matrix Done!"






