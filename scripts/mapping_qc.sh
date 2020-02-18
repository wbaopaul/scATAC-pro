#!/bin/bash

## output mapping qc results

input_dir=$1  ## the directory where the mapping results are

curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

output_dir=${OUTPUT_DIR}/summary
mkdir -p $output_dir

ncore=$(nproc --all)
ncore=$(($ncore - 1))

input_pre=${input_dir}/${OUTPUT_PREFIX}
output_pre=${output_dir}/${OUTPUT_PREFIX}

${SAMTOOLS_PATH}/samtools flagstat -@ $ncore ${input_pre}.positionsort.bam > ${output_pre}.flagstat.txt
${SAMTOOLS_PATH}/samtools idxstats -@ $ncore ${input_pre}.positionsort.bam > ${output_pre}.idxstat.txt

if [[ ! -f ${input_pre}.positionsort.MAPQ${MAPQ}.bam ]];then
	${SAMTOOLS_PATH}/samtools view -@ $ncore -h -b -q ${MAPQ} ${input_pre}.positionsort.bam > ${input_pre}.positionsort.MAPQ${MAPQ}.bam
fi

if [[ ! -f ${input_pre}.positionsort.MAPQ${MAPQ}.bam.bai ]];then
	${SAMTOOLS_PATH}/samtools index -@ $ncore ${input_pre}.positionsort.MAPQ${MAPQ}.bam
fi

${SAMTOOLS_PATH}/samtools flagstat -@ $ncore ${input_pre}.positionsort.MAPQ${MAPQ}.bam > ${output_pre}.MAPQ${MAPQ}.flagstat.txt
${SAMTOOLS_PATH}/samtools idxstats -@ $ncore ${input_pre}.positionsort.MAPQ${MAPQ}.bam > ${output_pre}.MAPQ${MAPQ}.idxstat.txt


tmp_sam_file=${output_dir}/tmp.sam
${SAMTOOLS_PATH}/samtools view -@ $ncore -q 5 -f 0x2 ${input_pre}.positionsort.bam > $tmp_sam_file

if [ $MAPPING_METHOD == 'bwa' ]; then
   total_uniq_mapped=$( wc -l ${tmp_sam_file} | cut -d ' ' -f1 )  ## number of unique mapped reads
elif [ $MAPPING_METHOD == 'bowtie' ]; then
   total_uniq_mapped=$( grep -E "@|NM:" $tmp_sam_file | grep -v "XS:" | wc -l )
else 
   total_uniq_mapped=$( grep -E "@|NM:" $tmp_sam_file | grep -v "XS:" | wc -l )
fi


total_uniq_mapped=$((${total_uniq_mapped}/2))

total_pairs=$(grep 'paired in' ${output_pre}.flagstat.txt | cut -d ' ' -f1) 
total_pairs=$((${total_pairs}/2)) 
total_pairs_mapped=$(grep 'with itself and mate mapped'  ${output_pre}.flagstat.txt | cut -d ' ' -f1)
total_pairs_mapped=$((${total_pairs_mapped}/2)) 
total_mito_mapped=$(grep chrM ${output_pre}.idxstat.txt | cut -f3)
total_mito_unmapped=$(grep chrM ${output_pre}.idxstat.txt | cut -f4)
total_mito=$((${total_mito_mapped}/2 + ${total_mito_unmapped}/2))
total_mito_mapped=$((${total_mito_mapped}/2))
total_dups=$(grep 'duplicates' ${output_pre}.flagstat.txt | cut -d ' ' -f1)
total_dups=$((${total_dups}/2))



total_pairs_MAPQH=$(grep 'with itself and mate mapped'  ${output_pre}.MAPQ${MAPQ}.flagstat.txt | cut -d ' ' -f1)
total_pairs_MAPQH=$((${total_pairs_MAPQH}/2)) 
total_mito_MAPQH=$(grep chrM ${output_pre}.MAPQ${MAPQ}.idxstat.txt | cut -f3)
total_mito_MAPQH=$((${total_mito_MAPQH}/2))
total_dups_MAPQH=$(grep 'duplicates' ${output_pre}.MAPQ${MAPQ}.flagstat.txt | cut -d ' ' -f1)
total_dups_MAPQH=$((${total_dups_MAPQH}/2))

rm ${output_pre}.idxstat.txt 
rm ${output_pre}.flagstat.txt 

rm ${output_pre}.MAPQ${MAPQ}.idxstat.txt 
rm ${output_pre}.MAPQ${MAPQ}.flagstat.txt 

#print to file
echo "Total_Pairs    $total_pairs" > ${output_pre}.MappingStats 
echo "Total_Pairs_Mapped    $total_pairs_mapped" >> ${output_pre}.MappingStats 
echo "Total_Uniq_Mapped    $total_uniq_mapped" >> ${output_pre}.MappingStats 
#echo "Total_Mito    $total_mito" >> ${output_pre}.MappingStats 
echo "Total_Mito_Mapped    $total_mito_mapped" >> ${output_pre}.MappingStats 
echo "Total_Dups    $total_dups" >> ${output_pre}.MappingStats 


echo "Total_Pairs_MAPQ${MAPQ}    $total_pairs_MAPQH" >> ${output_pre}.MappingStats 
echo "Total_Mito_MAPQ${MAPQ}    $total_mito_MAPQH" >> ${output_pre}.MappingStats 
echo "Total_Dups_MAPQ${MAPQ}    $total_dups_MAPQH" >> ${output_pre}.MappingStats 

rm $tmp_sam_file



