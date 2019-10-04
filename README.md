scATAC-pro
================

A comprehensive pipeline for single cell ATAC-seq data processing and analysis

Workflow
--------

scATAC-pro incorporates two main steps, preprocessing and downstream analysis. The preprocessing step takes raw fastq files as input and outputs peak-by-cell count matrix. It consists of demultiplexing, adaptor trimming, mapping, peak calling, cell calling, signal generating and quality controlling modules. The downstream analysis is comprised of dimension reduction, cell clustering, differential accessibility analysis, TF motif enrichment analysis and footprinting analysis. We provide flexible options for most of the modules.

<center>
<img src="/Users/wenbaoyu/fig1_v2.png"
     alt="Workflow"
     style="margin-right: 5px; height: 15cm; width: 15cm" />
</center>


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

### Tools users should install

-   R (&gt;=3.6.0)
-   Python (&gt;=3.6.0)
-   bowtie/bowtie2 (if users don't want to use bwa aligner)

### Tools required

**Will be automatically installed if NOT detected.**

-   BWA (&gt;=0.7.17)
-   MACS2 (&gt;=2.1.1)
-   samtools (&gt;=1.9)
-   bedtools (&gt;=2.27.1)
-   deepTools (&gt;=3.2.1)
-   trim\_galore (&gt;=0.6.3)
-   R packaages: devtools, flexdashboard, png, data.table, Matirx, Rcpp, ggplot2, flexmix, optparse, magrittr, readr, Seurat, bedr, gridExtra, ggrepel, kableExtra, viridis, RColorBrewer,pheatmap,motifmatchr, chromVAR, chromVARmotifs, SummarizedExperiment, BiocParallel

### Tools for additional modules or options

-   RGT (for footprint analysis, will ask whether you want to install it since the installation is done through conda, which takes a while and you may not want to conduct footprint analysis)
-   Trimmomatic
-   R packages (DESeq2, cisTopic, RcisTarget, AUCell, BSgenome.Hsapiens.UCSC.hg38, BSgenome.Hsapiens.UCSC.hg19, BSgenome.Mmusculus.UCSC.mm10, BSgenome.Mmusculus.UCSC.mm9, clusterProfiler, VisCello)

Quick start
-----------

-   **IMPORTANT**: the paramters and options should be specified in a configure file in text format: Copy and edit the configure\_user.txt file in this repository and then in your terminal run:

<!-- -->

    $ scATAC-pro -s process 
                 -i pe1_fastq,pe2_fastq,index_fastq 
                 -c configure_user.txt 

    $ scATAC-pro -s downstream 
                 -i output/filtered_matrix/YOUR_CELL_CALLER/matrix.mtx 
                 -c configure_user.txt
    ## YOUR_CELL_CALLER is specified in your configure_user.txt file

-   For data processing, if fastq files have been demultipled as the required format: the barcode was recorded in the name of each read like @barcode:ORIGIN\_READ\_NAME , you can skip the demultipling step by running

<!-- -->

    $ scATAC-pro -s process_no_dex 
                 -i pe1_fastq,pe2_fastq
                 -c configure_user.txt 

-   The **output** will be saved under ./output as default

Run scATAC-pro step by step
---------------------------

-   **IMPORTANT**: you can run the pipeline sequentially. The input of the later step was saved in output of previous step. The following commands using fastq files downloaded from [PBMC10k 10X Genomics](https://support.10xgenomics.com/single-cell-atac/datasets/1.1.0/atac_v1_pbmc_10k?) as an example:

-   *Combine data from different lanes*

<!-- -->


    $ cat atac_pbmc_10k_v1_S1_L001_R1_001.fastq.gz atac_pbmc_10k_v1_S1_L002_R1_001.fastq.gz > pe1_fastq

    $ cat atac_pbmc_10k_v1_S1_L001_R3_001.fastq.gz atac_pbmc_10k_v1_S1_L002_R3_001.fastq.gz > pe2_fastq

    $ cat atac_pbmc_10k_v1_S1_L001_R2_001.fastq.gz atac_pbmc_10k_v1_S1_L002_R2_001.fastq.gz > index_fastq

-   *Run the pipeline sequentially*

<!-- -->

    $ scATAC-pro -s demplx_fastq 
                 -i pe1_fastq,pe2_fastq,index_fastq 
                 -c configure_user.txt 

    $ scATAC-pro -s trimming 
                 -i output/demplxed_fastq/demplxed_pe1_fastq,
                    output/demplxed_fastq/demplxed_pe2_fastq
                 -c configure_user.txt 


    $ scATAC-pro -s mapping 
                  -i output/trimmed_fastq/trimmed_pe1_fastq,
                     output/trimmed_fastq/trimmed_pe2_fastq 
                  -c configure.txt 

    $ scATAC-pro -s call_peak 
                 -i output/mapping_result/pbmc10k.positionsort.MAPQ30.bam
                 -c configure.txt 

    $ scATAC-pro -s aggr_signal 
                 -i output/mapping_result/pbmc10k.positionsort.MAPQ30.bam 
                 -c configure.txt 
                 
    $ scATAC-pro -s get_mtx 
                 -i output/summary/pbmc10k.fragments.bed 
                 -c configure.txt 

    $ scATAC-pro -s qc_per_barcode 
                 -i output/summary/pbmc10k.fragments.bed 
                 -c configure.txt

    $ scATAC-pro -s call_cell
                 -i output/raw_matrix/YOUR_PEAK_CALLER/matrix.mtx
                 -c configure.txt
                 
    $ scATAC-pro -s clustering
                 -i output/filtered_matrix/YOUR_CELL_CALLER/matrix.mtx 
                 -c configure_user.txt

    $ scATAC-pro -s motif_analysis
                 -i output/filtered_matrix/YOUR_CELL_CALLER/matrix.mtx 
                 -c configure_user.txt
                 
    $ scATAC-pro -s split_bam
                 -i output/downstream_analysis/YOUR_CELL_CALLER/cell_cluster_table.txt
                 -c configure_user.txt

    $ scATAC-pro -s footprint
                 -i output/downstream_analysis/YOUR_CELL_CALLER/data_by_cluster/cluter_0.bam,
                    output/downstream_analysis/YOUR_CELL_CALLER/data_by_cluster/cluter_1.bam
                 -c configure_user.txt
                 
    $ scATAC-pro -s runDA
                 -i output/filtered_matrix/YOUR_CELL_CALLER/seurat_obj.rds 
                 -c configure_user.txt
                 
    $ scATAC-pro -s runGO
                 -i output/filtered_matrix/YOUR_CELL_CALLER/differential_peak_cluster_table.txt 
                 -c configure_user.txt
      
                 
    $ scATAC-pro -s report
                 -i output/summary
                 -c configure_user.txt
                 
    ## perform integrated analysis             
    $ scATAC-pro -s integrate
                 -i bam_file1,bam_file2,(bam_file3...)
                 -c configure_user.txt
                

Detailed Usage
--------------

    $ scATAC-pro --help
    usage : scATAC-pro -s STEP -i INPUT -c CONFIG [-o] [-h] [-v]
    Use option -h|--help for more information

    scATAC-pro 1.0.0
    ---------------
    OPTIONS

       [-s|--step ANALYSIS_STEP] : run a analytic step (or combinatorial steps) of the scATAC-pro workflow, supportting steps:
          dex_fastq: perform demultiplexing
                               input: fastq files for both reads and index, separated by comma like:
                                      fastq1,fastq2,index_fastq1,index_fastq2, index_fastq3...;
                               output: demultiplexed fastq1 and fastq2 files 
          mapping: perform reads alignment
                             input: fastq files, separated by comma for each paired end
                             output: position sorted bam file, mapping qc stat and fragment.bed
          call_peak: call peaks for aggregated data
                               input: BAM file path
                               output: peaks in bed file format
          get_mtx: build raw peak by barcode matrix
                             input: BAM file path
                             output: fragment.bed file and sparse matrix in Matrix Market format
          aggr_signal: generate aggregated signal, which can be upload to and view
                                 in genome browser
                                 input: require BAM file path
                                 output: bw and bedgraph file
          qc_per_barcode: quality control per barcode
                                    input: fragment.bed file
                                    output: qc_per_barcode.summary
          process: processing data - including dex_fastq, mapping, call_peak, get_mtx,
                                aggr_signal, qc_per_barcode abd call_cell
                                input: fastq files for both reads and index, separated by comma like:
                                       fastq1,fastq2,index_fastq1,index_fastq2, index_fastq3...; 
                                output: cell peak matrix and all intermediate results 
          call_cell: cell calling
                               input: raw peak barcode sparse matrix file path
                               output: filtered peak by cell matrix
          clustering: cell clustering
                               input: filtered peak by cell matrix file path
                               output: seurat objects with clustering label in the metadata (.rds file) and 
                                       barcodes with cluster labels (bed file)
          motif_analysis: doing motif analysis
                               input: filtered peak by cell matrix file path
                               output: TF by cell matrix indicating TF accessibility (chromVAR object)
          split_bam: split bam file into different clusters
                               input: barcodes with cluster label (.txt file, outputed from clustering)
                               output: .bam (saved under downstream/CELL_CALLER/data_by_cluster), .bw, .bedgr                                                            aph (save under output/signal/) file for each cluster
          footprint: doing footprinting analysis
                               input: bam files of two clusters, separated by comma like, bam1,bam2
                               output: footprint summary statistics (saved under output/downstream/CELL_CALLE                                                            R/footprinting/)
          integrate: doing integration analysis for two ore more samples
                               input: bam files, separated by comma like, bam1,bam2
                               output: save all intemediate results under output/integrated/)
          downstream: do all downstream analysis, including clustering, motif_analysis, 
                                split_bam (optional) and footprinting analysis (optional)
                                input: filtered matrix file path
                                output: all outputs from each step
          report: generate report in html file
                            input: directory to output report
                            output: summary report in html format
       -i|--input INPUT : input data, different types of input data are required for different steps;
       -c|--conf CONFIG : configuration file for parameters (if exists) for each step
       [-o|--output_dir : folder to save results, default output/ under the curret directory; sub-folder will be created automatically for each step
       [-h|--help]: help
       [-v|--version]: version
