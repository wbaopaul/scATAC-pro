#!/bin/bash

## call peaks on given cell barcodes 

input_bcs=$1


# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf $2
read_conf $3

out_prefix=${OUTPUT_PREFIX}.${MAPPING_METHOD}

input_bam0=${OUTPUT_DIR}/mapping_result/${out_prefix}.positionsort.MAPQ${MAPQ}.bam

## select bam given barcodes
${PERL_PATH}/perl ${curr_dir}/src/extract_bam4bcs.pl --barcode_file $inpput_bcs --bam_file $input_bam0 --output_prefix ${OUTPUT_DIR}/mapping_result/reads4selectedBarcodes --samtools_path $SAMTOOLS_PATH


input_bam=${OUTPUT_DIR}/mapping_result/reads4selectedBarcodes.bam


peaks_dir="${OUTPUT_DIR}/peaks"
mkdir -p $peaks_dir

output_prefix=${outptu_prefix}.selectedBarcodes


## call peaks
if [ ${PEAK_CALLER} = 'macs2' ];then
	echo "--Using macs2... "
	unset PYTHONPATH
	work_dir=${peaks_dir}/macs2
	mkdir -p $work_dir
    echo $out_prefix
	${MACS2_PATH}/macs2 callpeak -t $input_bam --outdir $work_dir -n $output_prefix -f BAM $MACS2_OPTS 

	## remove peaks overlapped with blacklist
	${BEDTOOLS_PATH}/bedtools intersect -a ${work_dir}/${out_prefix}_peaks.narrowPeak -b $BLACKLIST -v \
	    > ${work_dir}/${out_prefix}_peaks_BlacklistRemoved.bed
fi


if [ ${PEAK_CALLER} = 'GEM' ];then
	echo "--Using GEM..."
	work_dir=${peaks_dir}/GEM
	mkdir -p $work_dir
 	${SAMTOOLS_PATH}/samtools view -@ 3 -h $input_bam > $work_dir/tmp.sam
 	java -Xmx40G -jar ${GEM_PATH}/gem.jar --t 3 --g ${GEM_PATH}/hg38.chrom.sizes --s 2000000000 --d ${GEM_PATH}/Read_Distribution_default.txt --expt $input_bam --f SAM --out $work_dir/${out_prefix}
	sed -1 '1d' $work_dir/${out_prefix}/${out_prefix}.GPS_events.txt 
	# extend the binding site +/-150bp as peak
	cut -f1 ${wrok_dir}/${out_prefix}/${out_prefix}.GPS_events.txt | awk -F ":" '{print $1,"\t", $2}' | awk '{print "chr"$1, "\t", $2-150, "\t", $2+150}' > ${work_dir}/${out_prefix}_peaks.bed

	## remove peaks overlapped with blacklist
	${BEDTOOLS_PATH}/bedtools intersect -a ${work_dir}/${out_prefix}_peaks.bed -b $BLACKLIST -v \
	    > ${peaks_dir}/${out_prefix}_peaks_BlacklistRemoved.bed
fi

echo "Call peaks done !"
