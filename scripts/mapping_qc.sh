#!/bin/bash

## output mapping qc results

input_dir=$1  ## the directory where the mapping results are
output_dir=${2}/qc_result
mkdir -p $output_dir

ncore=$(nproc --all)
ncore=$(($ncore - 1))

input_pre=${input_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}
output_pre=${output_dir}/${OUTPUT_PREFIX}.${MAPPING_METHOD}

${SAMTOOLS_PATH}/samtools flagstat -@ $ncore ${input_pre}.positionsort.bam > ${output_pre}.flagstat.txt
${SAMTOOLS_PATH}/samtools idxstats -@ $ncore ${input_pre}.positionsort.bam > ${output_pre}.idxstat.txt

if [[ ! -f ${input_pre}.positionsort.MAPQ30.bam ]];then
	${SAMTOOLS_PATH}/samtools view -@ $ncore -h -b -q 30 ${input_pre}.positionsort.bam > ${input_pre}.positionsort.MAPQ30.bam
fi

if [[ ! -f ${input_pre}.positionsort.MAPQ30.bam.bai ]];then
	${SAMTOOLS_PATH}/samtools index -@ $ncore ${input_pre}.positionsort.MAPQ30.bam
fi

${SAMTOOLS_PATH}/samtools flagstat -@ $ncore ${input_pre}.positionsort.MAPQ30.bam > ${output_pre}.MAPQ30.flagstat.txt
${SAMTOOLS_PATH}/samtools idxstats -@ $ncore ${input_pre}.positionsort.MAPQ30.bam > ${output_pre}.MAPQ30.idxstat.txt


tmp_sam_file=${output_dir}/tmp.sam
${SAMTOOLS_PATH}/samtools view -@ $ncore -q 1 ${input_pre}.positionsort.bam > $tmp_sam_file

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



total_pairs_MAPQ30=$(grep 'with itself and mate mapped'  ${output_pre}.MAPQ30.flagstat.txt | cut -d ' ' -f1)
total_pairs_MAPQ30=$((${total_pairs_MAPQ30}/2)) 
total_mito_MAPQ30=$(grep chrM ${output_pre}.MAPQ30.idxstat.txt | cut -f3)
total_mito_MAPQ30=$((${total_mito_MAPQ30}/2))
total_dups_MAPQ30=$(grep 'duplicates' ${output_pre}.MAPQ30.flagstat.txt | cut -d ' ' -f1)
total_dups_MAPQ30=$((${total_dups_MAPQ30}/2))



#print to file
echo "total_pairs    $total_pairs" > ${output_pre}.summary 
echo "total_uniq_mapped    $total_uniq_mapped" >> ${output_pre}.summary 
echo "total_pairs_mapped    $total_pairs_mapped" >> ${output_pre}.summary 
echo "total_mito    $total_mito" >> ${output_pre}.summary 
echo "total_mito_mapped    $total_mito_mapped" >> ${output_pre}.summary 
echo "total_dups    $total_dups" >> ${output_pre}.summary 


echo "total_pairs_MAPQ30    $total_pairs_MAPQ30" >> ${output_pre}.summary 
echo "total_mito_MAPQ30    $total_mito_MAPQ30" >> ${output_pre}.summary 
echo "total_dups_MAPQ30    $total_dups_MAPQ30" >> ${output_pre}.summary 

rm $tmp_sam_file

## using R to make the result as a table/plot


echo "MAPPING QC Done!"
