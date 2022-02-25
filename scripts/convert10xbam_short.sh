#!/bin/bash

set -e

input_files=$1 ## specify input and output bam file names
input_files=(${input_files//,/ })

input_bam=${input_files[0]}
output_bam=${input_files[1]}
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

ncore=$(nproc --all)
ncore=$((${ncore} - 4))
if [[ $ncore -gt 30 ]]; then
    ncore=30
elif [[ $ncore -lt 1 ]]; then
    ncore=1
fi

# create a bam file with the barcode embedded into the read name
cat <( cat ${mapRes_dir}/${bamName}.header ) \
<( samtools view -@ $ncore $input_bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CB"], $0 }' ) \
| samtools view -@ $ncore -bS - > $output_bam

rm ${mapRes_dir}/${bamName}.header

${SAMTOOLS_PATH}/samtools index -@ $ncore $output_bam

