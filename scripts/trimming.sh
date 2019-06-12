#!/bin/bash



set -e

ff=$1
output_dir=$2
output_dir=${output_dir}/trimmed_fastq
mkdir -p $output_dir

fastqs=(${ff//,/ }) ## suppose the first fastq is the read file, the others are index fastq files
nfile=${#fastqs[@]}
kk=$(( $nfile ))

if [[ $kk < 2 ]];then
  echo -e "Erro: Provide at least two fastq files: with the first one as read fastq, the others as index fastq files" >&2
  exit
fi

curr_dir=`dirname $0`

ncore=$(nproc --all)
ncore=$(($ncore - 1))

## the first barcode was add to the read name after @, and : was used to concatenate to the original name
prefix0=${fastqs[0]}
prefix0=${prefix0##*/}

prefix1=$(basename ${fastqs[1]})
#prefix1=${fastqs[1]}
#prefix1=${prefix1##*/}

#java -jar ${Trimmomatic_PATH}/*jar PE ${fastqs[0]} ${fastqs[1]} ${output_dir}/trimmed_paired_${prefix0} ${output_dir}/trimmed_unpaired_${prefix0} \
#    ${output_dir}/trimmed_paired_${prefix1} ${output_dir}/trimmed_unpaired_${prefix1} \
#    ILLUMINACLIP:${Trimmomatic_PATH}/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:25

/mnt/isilon/cbmi/tan_lab/yuw1/local_tools/bin/trim_galore -j $ncore -o $output_dir  ${fastqs[0]} ${fastqs[1]}   --paired
echo "Trimming Done!"
