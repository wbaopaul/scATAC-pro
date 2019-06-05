#!/bin/bash

#### only keep proper paired, high mapq reads
#### for down stream analysis
#### with further filter barcodes with reads less than 10

input_bam=$1

filterBam_dir="${2}/filtered_bam"
mkdir -p $filterBam_dir

ncore=$(nproc --all)
ncore=$(($ncore - 1))

out_prefix=${OUTPUT_PREFIX}.${MAPPING_METHOD}

echo "filter by MAPQ${MAPQ}..."
${SAMTOOLS_PATH}/samtools view -@ $ncore -h -b -q $MAPQ $input_bam -o ${filterBam_dir}/tmp.bam

echo "only keep proper paired pairs..."
${SAMTOOLS_PATH}/samtools view -@ $ncore -f 0x2 -h -b  ${filterBam_dir}/tmp.bam -o ${filterBam_dir}/${out_prefix}.paired.MAPQ${MAPQ}.bam

rm ${filterBam_dir}/tmp.bam

echo "remove duplicates..."
${SAMTOOLS_PATH}/samtools markdup -@ $ncore -r ${filterBam_dir}/${out_prefix}.paired.MAPQ${MAPQ}.bam ${filterBam_dir}/${out_prefix}.dedup.paired.MAPQ${MAPQ}.bam 

rm ${filterBam_dir}/${out_prefix}.paired.MAPQ${MAPQ}.bam


${SAMTOOLS_PATH}/samtools view -h -@ $ncore  ${filterBam_dir}/${out_prefix}.dedup.paired.MAPQ${MAPQ}.bam > ${filterBam_dir}/tmp.sam



echo "further remove barcodes with less than 10 reads..."

curr_dir=`dirname $0`
perl ${curr_dir}/cal_frac_mito.pl --read_file ${filterBam_dir}/tmp.sam  --read_length 50  --output_file ${filterBam_dir}/tmp.bed
sed -i '1d' ${filterBam_dir}/tmp.bed
awk '$2 > 10 {print $1}' ${filterBam_dir}/tmp.bed > ${filterBam_dir}/barcodes.reads.GT10

perl ${curr_dir}/extract_sam_gbarcodes.pl --barcode_file ${filterBam_dir}/barcodes.reads.GT10 --read_file ${filterBam_dir}/tmp.sam --output_file ${filterBam_dir}/${out_prefix}.dedup.paired.MAPQ${MAPQ}.bcReadsGT10.sam

rm ${filterBam_dir}/tmp.sam
rm ${filterBam_dir}/tmp.bed


echo "Filtering Bam Done!"





