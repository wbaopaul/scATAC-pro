#!/bin/bash

input_bam=$1

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf $2
read_conf $3

peaks_dir="${OUTPUT_DIR}/peaks"
mkdir -p $peaks_dir

out_prefix=${OUTPUT_PREFIX}

## call peaks
if [ "${PEAK_CALLER}" = 'MACS2' ];then
	echo "--Using MACS2... "
	unset PYTHONHOME
	unset PYTHONPATH
	work_dir=${peaks_dir}/MACS2
	mkdir -p $work_dir
	${MACS2_PATH}/macs2 callpeak -t $input_bam --outdir $work_dir -n $out_prefix -f BAM $MACS2_OPTS 
	#${MACS2_PATH}/macs2 callpeak -t $input_bam --outdir $peaks_dir -f BAM $MACS2_OPTS --nomodel --extsize 147

	## remove peaks overlapped with blacklist
	${BEDTOOLS_PATH}/bedtools intersect -a ${work_dir}/${out_prefix}_peaks.narrowPeak -b $BLACKLIST -v \
	    > ${work_dir}/${out_prefix}_features_BlacklistRemoved.bed
fi


if [ "${PEAK_CALLER}" = 'MUSIC' ];then
	echo "--Using MUSIC..."
	work_dir=${peaks_dir}/MUSIC
	mkdir -p $work_dir/chip/sorted
	mkdir -p $work_dir/chip/dedup
	${SAMTOOLS_PATH}/samtools view $input_bam | ${MUSIC_PATH}/MUSIC -preprocess SAM stdin ${work_dir}/chip

	${MUSIC_PATH}/MUSIC -sort_reads ${work_dir}/chip ${work_dir}/chip/sorted
	${MUSIC_PATH}/MUSIC -remove_duplicates ${work_dir}/chip/sorted 2 ${work_dir}/chip/dedup
	${MUSIC_PATH}/MUSIC -get_multiscale_broad_ERs -chip ${work_dir}/chip/dedup \
        -mapp /mnt/isilon/tan_lab/yuw1/run_scATAC-pro/Mappability_MAP/map50_hg38 \
        -l_mapp 50 -begin_l 1000 -end_l 16000 -step 1.5
 

	## remove peaks overlapped with blacklist
fi



if [ "${PEAK_CALLER}" = 'GEM' ];then
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
	    > ${peaks_dir}/${out_prefix}_features_BlacklistRemoved.bed
fi

if [ "${PEAK_CALLER}" = 'COMBINED' ];then
    curr_dir=`dirname $0`
    ${curr_dir}/iter_peak.sh $1 $2 $3
fi

if [ "${PEAK_CALLER}" = 'BINNING' ];then
    echo "--Binning genome"
	work_dir=${peaks_dir}/BIN
	mkdir -p $work_dir
    bin_file=${mat_dir}/${output_prefix}_bin.bed
    ${BEDTOOLS_PATH}/bedtools makewindows -g $CHROM_SIZE_FILE -w $BIN_RESL > $bin_file
	${BEDTOOLS_PATH}/bedtools intersect -a ${work_dir}/${out_prefix}_bin.bed -b $BLACKLIST -v \
	    > ${peaks_dir}/${out_prefix}_features_BlacklistRemoved.bed
fi
   
echo "Call peaks done !"
