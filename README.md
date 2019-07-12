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

-   Python (&gt;=3.6.0)
-   R (&gt;=3.5.3)

### Tools will be installed during installation (if NOT detected)

-   BWA will be installed if none of BWA/bowtie/bowtie2 was detected
-   MACS2 (&gt;=2.1.1)
-   samtools (&gt;=1.9)
-   deepTools (&gt;=3.2.1)
-   trim\_galore (&gt;=0.6.3)
-   R packaages: data.table, ggplot2, Seurat, chromVAR, DropletUtils, Rcpp, flexmix, optparse, magrittr, readr, SummarizedExperiment, BiocParallel, motifmatchr, JASPAR2016
-   Python module: numpy

### Optional tools/packages for optional downstream analysis (will try to install it if USED and NOT detected)

-   Trimmomatic (if Trimmomatic was chosen rather than trim\_galore)
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

        [-s|--step ANALYSIS_STEP] : run a step (or combinatorial steps) of the scATAC-pro workflow, supportting steps:
          dex_fastq: perform demultiplexing
                     input: fastq files for both reads and index, separated by comma like:
                            fastq1,fastq2,index_fastq1,index_fastq2, index_fastq3...;
                     output: demultiplexed fastq1 and fastq2 files 
          mapping: perform reads alignment
                     input: fastq files, separated by comma for each paired end
                     output: position sorted bam file and mapping qc stat
          call_peak: call peaks for aggregated data
                     input: BAM file
                     output: peaks in bed file format
          get_mtx: build raw peak by barcode matrix
                     input: BAM file
                     output: fragment.bed file and sparse matrix in Matrix Market format
          aggr_signal: generate aggregated signal, which can be upload to and view
                       in genome browser
                     input: require BAM files
                     output: bw and bedgraph file
          qc_per_barcode: quality control per barcode
                     input: fragment.bed file
                     output: qc_per_barcode.summary
          preprocess: perform preprocessing - including dex_fastq, mapping, call_peak, get_mtx,
                      aggr_signal, qc_per_barcode abd call_cell
                     input: fastq files for both reads and index, separated by comma like:
                            fastq1,fastq2,index_fastq1,index_fastq2, index_fastq3...; 
                     output: cell peak matrix and all intermediate results 
          call_cell: cell calling
                     input: raw peak barcode sparse matrix file
                     output: filtered peak by cell matrix
          clustering: cell clustering
                     input: filtered peak by cell matrix file
                     output: seurat objects with clustering label in the metadata (.rds file) and 
                                       barcodes with cluster labels (bed file)
          motif_analysis: doing motif analysis
                     input: filtered peak by cell matrix file
                     output: TF by cell matrix indicating TF accessibility (chromVAR object)
          split_bam: split bam file into different clusters
                     input: barcodes with cluster label (bed file)
                     output: footprint summary statistics (hint output)
          footprint: doing footprinting analysis
                     input: bam files of two clusters, separated by comma like, bam1,bam2
                     output: footprint summary statistics (hint output)
          downstream: do all downstream analysis, including clustering, motif_analysis, 
                      split_bam (optional) and footprinting analysis (optional)
                     input: sparse matrix file
                     output: seurat object (in .rds format), differential TF accessibility score                           matrix and footprint summary
          report: generate report in html file
                     input: directory to output report
                     output: summary report in html format
       -i|--input INPUT : input data, different types of input data are required for different steps;
       -c|--conf CONFIG : configuration file for parameters (if exists) for each step
       [-o|--output_dir : folder to save results; sub-folder will be created automatically for each step
       [-p|--parallel] : if specified run scATAC-pro on a cluster
       [-h|--help]: help
       [-v|--version]: version

-   NOTE: the paramters and options are given in a text file: Copy and edit the 'configure\_user.txt' file
-   Detailed manual is here
-   Run scATAC-pro:

<!-- -->

    $ YOUR_INSTALL_PATH/scATAC-pro_1.0.0/scATAC-pro -s STEP -i INPUT -c configure_user.txt
