#!/bin/bash


#### run step by step #####


./scATAC-pro -s demplx_fastq -i test_fastq/test_R1.fastq,test_fastq/test_R2.fastq -c configure.txt -o ../tmp_output 
./scATAC-pro -s demplx_fastq -i test_fastq/test_R3.fastq,test_fastq/test_R2.fastq -c configure.txt -o ../tmp_output 
## note there should be no space after comma,


./scATAC-pro -s trimming -i test_fastq/test_R1.fastq,test_fastq/test_R3.fastq -c configure.txt -o ../tmp_output 

./scATAC-pro -s mapping -i ../tmp_output/demplxed_fastq/test_R1.fastq.dex,../tmp_output/demplxed_fastq/test_R3.fastq.dex -c configure.txt -o ../tmp_output > logs/map.log

./scATAC-pro -s filter_bam -i ../tmp_output/mapping_result/ss.markdup.bam -c configure.txt -o ../tmp_output > logs/filter_bam.log

./scATAC-pro -s call_peak -i ../tmp_output/filtered_bam/ss.dedup.MAPQ10.bam -c configure.txt -o ../tmp_output > logs/call_peak.log

./scATAC-pro -s get_peak_barcode_mat -i ../tmp_output/filtered_bam/ss.dedup.MAPQ10.bam -c configure.txt -o ../tmp_output > logs/get_peak_barcode_mat.log

./scATAC-pro -s generate_signal -i ../tmp_output/filtered_bam/ss.dedup.MAPQ10.bam -c configure.txt -o ../tmp_output > logs/generate_signal.log

./scATAC-pro -s qc_per_barcode -i ../tmp_output/filtered_bam/ss.dedup.MAPQ10.bam -c configure.txt -o ../tmp_output > logs/qc_per_barcode.log

