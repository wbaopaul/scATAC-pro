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

organism=hs
if [[ $GENOME_NAME =~ "mm" ]]; then
    organism=mm
fi

## call peaks
if [ "${PEAK_CALLER}" = 'MACS2' ];then
	echo "--Using MACS2... "
	unset PYTHONHOME
	unset PYTHONPATH
	work_dir=${peaks_dir}/MACS2
	mkdir -p $work_dir
	${MACS2_PATH}/macs2 callpeak -t $input_bam --outdir $work_dir -n $out_prefix -f BAM $MACS2_OPTS 
	#${MACS2_PATH}/macs2 callpeak -t $input_bam --outdir $peaks_dir -f BAM $MACS2_OPTS --nomodel --extsize 147

    ## remove peaks whose chromosome is not list in the chrom_size file
    ${R_PATH}/Rscript --vanilla ${curr_dir}/src/rmRedundantPeaks.R ${work_dir}/${out_prefix}_peaks.narrowPeak \
         $CHROM_SIZE_FILE
	
    ## remove peaks overlapped with blacklist
	${BEDTOOLS_PATH}/bedtools intersect -a ${work_dir}/${out_prefix}_peaks.narrowPeak -b $BLACKLIST -v \
	    > ${work_dir}/${out_prefix}_features_BlacklistRemoved.bed
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
    ${R_PATH}/Rscript --vanilla ${curr_dir}/src/rmRedundantPeaks.R ${work_dir}/${out_prefix}_peaks.bed $CHROM_SIZE_FILE
	${BEDTOOLS_PATH}/bedtools intersect -a ${work_dir}/${out_prefix}_peaks.bed -b $BLACKLIST -v \
	    > ${work_dir}/${out_prefix}_features_BlacklistRemoved.bed
fi

if [ "${PEAK_CALLER}" = 'COMBINED' ];then
    curr_dir=`dirname $0`
	work_dir=${peaks_dir}/COMBINED
    ${curr_dir}/iter_peak.sh $1 $2 $3
fi

if [ "${PEAK_CALLER}" = 'BIN' ];then
    echo "--Binning genome"
	work_dir=${peaks_dir}/BIN
	mkdir -p $work_dir
    bin_file=${work_dir}/${OUTPUT_PREFIX}_bin.bed
    ${BEDTOOLS_PATH}/bedtools makewindows -g $CHROM_SIZE_FILE -w $BIN_RESL > ${bin_file}
    ${R_PATH}/Rscript --vanilla ${curr_dir}/src/rmRedundantPeaks.R $bin_file $CHROM_SIZE_FILE
	${BEDTOOLS_PATH}/bedtools intersect -a ${bin_file} -b $BLACKLIST -v \
	    > ${work_dir}/${out_prefix}_features_BlacklistRemoved.bed
fi
  

## extend peaks which is shorter than 500bp 
${R_PATH}/Rscript --vanilla ${curr_dir}/src/extendPeaks.R ${work_dir}/${out_prefix}_features_BlacklistRemoved.bed
echo "Call peaks done !"
