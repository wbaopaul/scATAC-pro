#!/bin/bash

set -e
## input is the bam file and a txt file saved a cell barcode for each line, separated by comma
inputs=$1 
inputs=(${inputs//,/ })

input_bam=${inputs[0]}
ff=${inputs[1]}

# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

map_dir=${OUTPUT_DIR}/mapping_result
signal_dir=${OUTPUT_DIR}/signal
mkdir -p $map_dir
mkdir -p $signal_dir


curr_dir=`dirname $0`

## extract bam for cell and non-cells
${PERL_PATH}/perl ${curr_dir}/src/extract_bam4Cells.pl --cellbarcode_file $ff --bam_file $input_bam \
    --output_dir $map_dir --samtools_path $SAMTOOLS_PATH
${SAMTOOLS_PATH}/samtools view -@ 4 -bS ${map_dir}/cell_barcodes.sam > ${map_dir}/cell_barcodes.bam &
${SAMTOOLS_PATH}/samtools view -@ 4 -bS ${map_dir}/non_cell_barcodes.sam > ${map_dir}/non_cell_barcodes.bam &
echo "The bam file was split between cell and non-cell!"
wait

rm ${map_dir}/cell_barcodes.sam
rm ${map_dir}/non_cell_barcodes.sam

## QC using cell barcodes bam
CELL_MAP_QC=${CELL_MAP_QC^^}
if [[ "$CELL_MAP_QC"="TRUE" ]]; then

    echo "generate mapping stats for aggregated cell barcodes file..."
    bash ${curr_dir}/cell_mapping_qc.sh $map_dir $2 $3

    unset PYTHONHOME
    unset PYTHONPATH
    ncore=$(nproc --all)
    ncore=$((${ncore}/2))
    for file0 in $(find $map_dir -name "cell_barcodes.MAPQ${MAPQ}.bam")
    do
        echo "generate bw file..."
        pre=$(basename $file0)
        pre=${pre/.bam/}
        fname_bw=${signal_dir}/${pre}.bw
        ${SAMTOOLS_PATH}/samtools index -@ $ncore $file0
        ${DEEPTOOLS_PATH}/bamCoverage --numberOfProcessors max --normalizeUsing BPM \
             --bam $file0 --binSize 20 --skipNonCoveredRegions \
             --outFileName $fname_bw  
        echo "generate count around TSS..."
        ${DEEPTOOLS_PATH}/computeMatrix reference-point -S $fname_bw -R $TSS \
            -a 1200 -b 1200 -o ${signal_dir}/${pre}.aggregated.mtx.gz 
    done

    rm ${map_dir}/non_cell_barcodes.bam

    echo "Aggregated signal was generated for cell barcodes and non-cell barcodes!"
fi

