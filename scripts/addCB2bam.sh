#!/bin/bash

set -e

bam_file=$1

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

mapRes_dir=$( dirname $bam_file )
curr_dir=`dirname $0`

ncore=$(nproc --all)
ncore=$(($ncore - 1))

${SAMTOOLS_PATH}/samtools view -h -@ $ncore $bam_file |awk '{ if($0 ~ "^@") {print $0} else {split($1,a,":"); print $0"\tCB:Z:"a[1]} }' > ${mapRes_dir}/tmp.sam

new_file_name=${bam_file/.bam/_withCBtag.bam}
${SAMTOOLS_PATH}/samtools view -bS -@ $ncore ${mapRes_dir}/tmp.sam  > $new_file_name

rm ${mapRes_dir}/tmp.sam



