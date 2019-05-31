#!/bin/bash

## output mapping qc results

input_dir=$1
output_dir=$2
ncore=$(nproc --all)
ncore=$(($ncore - 1))

fname=${2}/${OUTPUT_PREFIX}.${MAPPING_METHOD}

${SAMTOOLS_PATH}/samtools flagstat -@ $ncore $input_dir/${fname}.markdup.bam > ${fname}.flagstat.txt
${SAMTOOLS_PATH}/samtools idxstats -@ $ncore $input_dir/${fname}.markdup.bam > ${fname}.idxstat.txt

${SAMTOOLS_PATH}/samtools flagstat -@ $ncore $input_dir/${fname}.MAPQ30.bam > ${fname}.MAPQ30.flagstat.txt
${SAMTOOLS_PATH}/samtools idxstats -@ $ncore $input_dir/${fname}.MAPQ30.bam > ${fname}.MAPQ30.idxstat.txt


tmp_sam=${2}/tmp.sam
${SAMTOOLS_PATH}/samtools view -@ $ncore -q 1 $input_dir/${fname}.markdup.bam > $tmp_sam

if [ $MAPPING_METHOD == 'bwa' ]; then
  total_uniq_map=$(wc -l $tmp_sam)  ## number of unique mapped reads
elif [ $MAPPING_METHOD == 'bowtie' ]; then
  total_uniq_map=$( grep -E "@|NM:" $tmp_sam | grep -v "XS:" | wc -l )
else 
  total_uniq_map=$( grep -E "@|NM:" $tmp_sam | grep -v "XS:" | wc -l )
fi

total_uniq_map=$(($total_uniq_map/2))
total_pairs=$(grep 'paired in' ${fname}.flagstat.txt | cut -d ' ' -f1) 
total_pairs=$((${total_pairs}/2)) 
total_pairs_mapped=$(grep 'properly paired'  ${fname}.flagstat.txt | cut -d ' ' -f1)
total_pairs_mapped=$((${total_pairs_mapped}/2)) 
total_mito_mapped=$(grep chrM ${fname}.idxstat.txt | cut -f3)
total_mito_unmapped=$(grep chrM ${fname}.idxstat.txt | cut -f4)
total_mito=$((${total_mito_mapped}/2 + ${total_mito_unmapped}/2))
total_dups=$(grep 'duplicates' ${fname}.flagstat.txt | cut -d ' ' -f1)
total_dups=$(($total_dups/2))


total_pairs30=$(grep 'paired in' ${fname}.MAPQ30.flagstat.txt | cut -d ' ' -f1) 
total_pairs30=$((${total_pairs30}/2)) 
total_pairs_mapped30=$(grep 'properly paired'  ${fname}.MAPQ30.flagstat.txt | cut -d ' ' -f1)
total_pairs_mapped30=$((${total_pairs_mapped30}/2)) 
total_mito_mapped30=$(grep chrM ${fname}.MAPQ30.idxstat.txt | cut -f3)
total_mito_mapped30=$((${total_mito_mapped30}/2))
total_dups30=$(grep 'duplicates' ${fname}.MAPQ30.flagstat.txt | cut -d ' ' -f1)
total_dups30=$((${total_dups30}/2))



#print to file
echo "total_uniq_map    $total_uniq_map" > ${fname}.summary 
echo "total_q30 $total_q30" >> ${fname}.summary 


rm $tmp_sam

## using R to make the result as a table/plot


