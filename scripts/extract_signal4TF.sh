#!/bin/bash

input_bw=$1 
input_tfbs=$2  ## a bed file (footprint result overlapping with motif)
tf_name=$3 ## a tf name, like CTCF

unset PYTHONPATH
DEEPTOOLS_PATH=$(which computeMatrix)
DEEPTOOLS_PATH=$(dirname $DEEPTOOLS_PATH)
echo "generate count around A TF..."

grep -w $tf_name $input_tfbs > tmp.bed
tfbs_dir=$(dirname $input_tfbs)
outname=${tfbs_dir}/${tf_name}_signal
${DEEPTOOLS_PATH}/computeMatrix reference-point -S $input_bw -R tmp.bed \
    -a 120 -b 120 -o ${outname}.mtx.gz


#average the signal and plot
