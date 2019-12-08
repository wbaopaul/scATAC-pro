#!/bin/bash

set -e

input_bam=$1

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

mapRes_dir="${OUTPUT_DIR}/mapping_result"
mkdir -p $mapRes_dir
curr_dir=`dirname $0`

echo "write the barcode in format of scATAC-pro:"
# extract the header file
bamName=$(basename $input_bam)
${SAMTOOLS_PATH}/samtools view $input_bam -H > ${mapRes_dir}/${bamName}.header

# create a bam file with the barcode embedded into the read name
cat <( cat ${mapRes_dir}/${bamName}.header ) \
<( samtools view $input_bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CB"], $0 }' ) \
| samtools view -bS - > ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.bam

rm ${mapRes_dir}/${bamName}.header


ncore=$(nproc --all)
ncore=$(($ncore - 1))

${SAMTOOLS_PATH}/samtools index -@ $ncore ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.bam 


## filtering low quality and/or deplicates for downstreame analysis
${SAMTOOLS_PATH}/samtools view -f 0x2 -b -h -q 30 -@ $ncore ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.bam -o ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.MAPQ30.bam 
${SAMTOOLS_PATH}/samtools index -@ $ncore ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.MAPQ30.bam 

if [ $MAPQ -ne 30 ]; then
     ${SAMTOOLS_PATH}/samtools view -f 0x2 -b -h -q $MAPQ -@ $ncore ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.bam -o ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.MAPQ${MAPQ}.bam 
fi


## mapping stats
echo "Summarizing mapping stats ..."

curr_dir=`dirname $0`
qc_dir=${OUTPUT_DIR}/summary
mkdir -p $qc_dir
bash ${curr_dir}/mapping_qc.sh ${mapRes_dir}  $2 $3



echo "Simple mapping stats summary Done!"

echo "Getting txt file for read pair (fragment) information"
${PERL_PATH}/perl ${curr_dir}/src/simply_bam2frags.pl --read_file ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.MAPQ${MAPQ}.bam \
        --output_file ${qc_dir}/${OUTPUT_PREFIX}.fragments.txt --samtools_path $SAMTOOLS_PATH

#echo "Remove duplicates"
#${SAMTOOLS_PATH}/samtools markdup -@ $ncore -r ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.MAPQ30.bam ${mapRes_dir}/${OUTPUT_PREFIX}.positionsort.MAPQ30.noDuplicates.bam 
#rm ${mapRes_dir}/${OUTPUT_PREFIX}.bam*
