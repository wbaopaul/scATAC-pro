#!/bin/bash



set -e

ff=$1
output_dir=$2
output_dir=${output_dir}/demplxed_fastq
mkdir -p $output_dir

fastqs=(${ff//,/ })
nfile=${#fastqs[@]}
kk=$(( $nfile-2 ))

curr_dir=`dirname $0`


for (( i==0; i<=$kk; i++ ))
do
  dex_prefix=${fastqs[$i]}
  dex_prefix=${dex_prefix##*/}
  ${PYTHON_PATH}/python ${curr_dir}/dex_fastq.py ${fastqs[$i]} ${output_dir}/demplxed_${dex_prefix}  ${fastqs[$nfile - 1]}
done

echo "Demultiplexing Done!"
