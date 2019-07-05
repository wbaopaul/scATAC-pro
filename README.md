scATAC-pro
================

A single cell ATAC-seq process pipeline

Introduction
------------

Installation
------------

-   Run the following command under your terminal, scATAC-pro will be installed in YOUR\_INSTALL\_PATH/scATAC-pro\_1.0.0

<!-- -->

    $ git clone https://github.com/wbaopaul/scATAC-pro.git
    $ cd scATAC-pro
    $ make configure prefix=YOUR_INSTALL_PATH
    $ make install
     

Dependencies
------------

### Tools user should install

-   Python3
-   R (&gt;=3.5.3)
-   MACS2 (&gt;=2.1.1)

### Tools will be installed during installation (if NOT detected)

-   BWA will be installed if none of BWA/bowtie/bowtie2 was detected
-   samtools (&gt;=1.9)
-   deepTools (&gt;=3.2.1)
-   trim\_galore (&gt;=0.6.3)

### R Rackages will be installed during installation (if NOT detected)

-   data.table, ggplot2, Seurat, chromVAR, DropletUtils, Rcpp, flexmix, optparse, magrittr, readr, SummarizedExperiment, BiocParallel, motifmatchr, JASPAR2016

### Optional tools/packages for optional downstream analysis (will try to install it if USED and NOT detected)

-   RGT (for footprint analysis)
-   R packages (cisTopic, RcisTarget, AUCell, BSgenome.Hsapiens.UCSC.hg38, BSgenome.Hsapiens.UCSC.hg19, BSgenome.Mmusculus.UCSC.mm10, BSgenome.Mmusculus.UCSC.mm9)

Usage
-----

    $ YOUR_INSTALL_PATH/scATAC-pro_1.0.0/scATAC-pro --help
    usage : YOUR_INSTALL_PATH/scATAC-pro_1.0.0/scATAC-pro -s STEP -i INPUT -c CONFIG [-o] [-p] [-h] [-v]
    Use option -h|--help for more information

    scATAC-pro 1.0.0
    ---------------
    OPTIONS

       [-s|--step ANALYSIS_STEP] : run only a subset of the scATAC-pro workflow; follow steps and cominatorial steps are provided:
          dex_fastq: perform demultiplexing - require fastq files both for reads and index, separated by comma like: fastq1,fastq2,index_fastq1,index_fastq2, index_fastq3...; output demultiplexed fastq1 and fastq2 files 
          mapping: perform reads alignment - require fastq files, separated by comma for each paired end
          call_peak: call peaks for aggregated data - require BAM files
          get_mtx: build raw peak by barcode matrix - require BAM files, output fragment.bed file and sparese matrix file
          aggr_signal: generate aggregated signal in bw file and bedgraph file, which can be upload and view in as genome browser - require BAM files, output bw and bedgraph file
          qc_per_barcode: quality control per barcode - require fragment.bed files
          call_cell: cell calling - require sparse matrix files
          preprocess: perform preprocessing - require fastq files both for reads and index, separated by comma like: fastq1,fastq2,index_fastq1,index_fastq2, index_fastq3...; output cell peak matrix and all intermediate results 
          downstream: do all downstream analysis, including clustering, motif_analysis, and footprinting analysis (optional) - require BAM files
          report: generate report in html file
       -i|--input INPUT : input data, different types of input data are required for different steps;
       -c|--conf CONFIG : configuration file for parameters (if exists) for each step
       [-o|--output_dir : folder to save results; sub-folder will be created automatically for each step
       [-p|--parallel] : if specified run scATAC-pro on a cluster
       [-h|--help]: help
       [-v|--version]: version

-   Copy and edit the 'configure\_user.txt' file
-   Detailed manual is here
-   Run scATAC-pro:

<!-- -->

    $ YOUR_INSTALL_PATH/scATAC-pro_1.0.0/scATAC-pro -s STEP -i INPUT -c configure_user.txt
