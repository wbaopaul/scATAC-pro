scATAC-pro
=================

A comprehensive pipeline for single cell ATAC-seq data processing and analysis


   * [scATAC-pro](#scatac-pro)
      * [Workflow](#workflow)
      * [Updates](#updates)
      * [Installation](#installation)
      * [Dependencies](#dependencies)
         * [Tools users should install](#tools-users-should-install)
         * [Tools required](#tools-required)
         * [Tools for additional modules or options](#tools-for-additional-modules-or-options)
      * [Quick start](#quick-start)
      * [Run scATAC-pro step by step](#run-scatac-pro-step-by-step)
      * [Detailed Usage](#detailed-usage)
      * [Run through docker or singularity](#run-through-docker-or-singularity)
      * [FAQs](#FAQs)
      * [Citation](#citation)


Workflow
--------

scATAC-pro incorporates two main steps, preprocessing and downstream analysis. The preprocessing step takes raw fastq files as input and outputs peak-by-cell count matrix. It consists of demultiplexing, adaptor trimming, mapping, peak calling, cell calling, signal generating and quality controlling modules. The downstream analysis is comprised of dimension reduction, cell clustering, differential accessibility analysis, TF motif enrichment analysis and footprinting analysis. We provide flexible options for most of the modules.


<p align="center">
  <img src="doc/fig1.png" width="480" title="">
</p>

Installation
------------

-   Note: you don't have to install it, you can use the docker or sigularity version if you prefer (see [Run through docker or singularity](#run-through-docker-or-singularity) ), note that the docker version was still under version 1.0.0 and will be updated to v1.1.0 soon
-   Run the following command under your terminal, scATAC-pro will be installed in YOUR\_INSTALL\_PATH/scATAC-pro\_1.1.0

<!-- -->

    $ git clone https://github.com/wbaopaul/scATAC-pro.git
    $ cd scATAC-pro
    $ make configure prefix=YOUR_INSTALL_PATH
    $ make install
     
Updates
------------
- current version: 1.1.0
- Feb, 2020
    * new module *visualize*, allowing interactively visualize and analysis the data
    * install rgt-hint (for footprint analysis) through miniconda3
    * add module *process_with_bam*, allowing process from aggragated bam file
    * integrate from peaks files, assume each sample was processed through scATAC-pro;
        output matrix with the same merged peaks/features and the previously called cells, along
        with a integrated seurat object
    * add new parameters in configure file: Top_Variable_Features, REDUCTION, nREDUCTION
    * enable all clustering methods mentioned in the manuscript, along with kmeans on PCs
    * file path changed to like downstreame_analysis/PEAK_CALLER/CELL_CALLER/..., indicating peak caller
    * qc_per_barcode requires too input files, separated by comma, see example and detailed usage
- Jan11, 2020 
    * add a new module mergePeaks to merge different peak files called from different samples or conditions
    * add a new module to reconstruct peak-cell matrix given a peak file, a fragment file and a barcodes.txt file
- Dec22, 2019 
    * corrected an error arised from using older version of chromVAR
- Dec11, 2019 
    * corrected a bug for demultiplexing multiple index files
- Dec7, 2019 
    * added a module convert10xbam to convert 10x position sorted bam file to scATAC-pro style
- Dec3, 2019 
    * updated module get_bam4Cells, with required inputs as a bam file and a txt file of barcodes, separated by comma







Dependencies
------------

### Tools users should install

-   R (&gt;=3.6.0)
-   Python (&gt;=3.6.0)

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
-   bowtie/bowtie2 (user install them if don't want to use bwa)
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
                  -c configure_user.txt 

    $ scATAC-pro -s call_peak 
                 -i output/mapping_result/pbmc10k.positionsort.MAPQ30.bam
                 -c configure_user.txt 

    $ scATAC-pro -s aggr_signal 
                 -i output/mapping_result/pbmc10k.positionsort.MAPQ30.bam 
                 -c configure_user.txt 
                 
    $ scATAC-pro -s get_mtx 
                 -i output/peaks/MACS2/pbmc10k_features_BlacklistRemoved.bed 
                 -c configure_user.txt 

    $ scATAC-pro -s qc_per_barcode 
                 -i output/summary/pbmc10k.fragments.txt,output/peaks/MACS2/pbmc10k_features_BlacklistRemoved.bed 
                 -c configure_user.txt

    $ scATAC-pro -s call_cell
                 -i output/raw_matrix/YOUR_PEAK_CALLER/matrix.mtx
                 -c configure_user.txt
                 
    $ scATAC-pro -s get_bam4Cells
                 -i output/mapping_result/pbmc10k.positionsort.bam,
                    output/filtered_matrix/YOUR_PEAK_CALLER/YOUR_CELL_CALLER/barcodes.txt
                 -c configure_user.txt

    $ scATAC-pro -s clustering
                 -i output/filtered_matrix/YOUR_PEAK_CALLER/YOUR_CELL_CALLER/matrix.mtx 
                 -c configure_user.txt

    $ scATAC-pro -s motif_analysis
                 -i output/filtered_matrix/YOUR_PEAK_CALLER/YOUR_CELL_CALLER/matrix.mtx 
                 -c configure_user.txt
                 
    $ scATAC-pro -s split_bam
                 -i output/downstream_analysis/YOUR_PEAK_CALLER/YOUR_CELL_CALLER/cell_cluster_table.txt
                 -c configure_user.txt

    $ scATAC-pro -s footprint
                 -i output/downstream_analysis/YOUR_PEAK_CALLER/YOUR_CELL_CALLER/data_by_cluster/cluter_0.bam,
                    output/downstream_analysis/YOUR_PEAK_CALLER/YOUR_CELL_CALLER/data_by_cluster/cluter_1.bam
                 -c configure_user.txt
                 
    $ scATAC-pro -s runDA
                 -i output/filtered_matrix/YOUR_PEAK_CALLER/YOUR_CELL_CALLER/seurat_obj.rds 
                 -c configure_user.txt
                 
    $ scATAC-pro -s runGO
                 -i output/filtered_matrix/YOUR_PEAK_CALLER/YOUR_CELL_CALLER/differential_peak_cluster_table.txt 
                 -c configure_user.txt
      
                 
    $ scATAC-pro -s report
                 -i output/summary
                 -c configure_user.txt
                 
    ## perform integrated analysis            
    $ scATAC-pro -s integrate
                 -i bam_file1,bam_file2,(bam_file3...)
                 -c configure_user.txt
                
    ## merge peaks that are within 200bp distance of each other            
    $ scATAC-pro -s mergePeaks
                 -i peak_file1,peak_file2,(peak_file3...),200
                 -c configure_user.txt

- After clustering, user can interactively visualize the cis-element chromatin accessibility using [VisCello](https://github.com/qinzhu/VisCello), in R:

```
scATAC-pro -s visualize -i output/downstream_analysis/YOUR_PEAK_CALLER/YOUR_CELL_CALLER/VisCello_obj -c configure_user.txt

```

Detailed Usage
--------------

    $ scATAC-pro --help
    usage : scATAC-pro -s STEP -i INPUT -c CONFIG [-o] [-h] [-v]
    Use option -h|--help for more information

    scATAC-pro 1.0.0
    ---------------
    OPTIONS

       [-s|--step ANALYSIS_STEP] : run an analytic step (or combinatorial steps) of the scATAC-pro workflow, supportting steps:
          demplx_fastq: perform demultiplexing
                               input: fastq files for both reads and index, separated by comma like:
                                      PE1_fastq,PE2_fastq,index1_fastq,inde2_fastq,index3_fastq...;
                                      differnet index will be embedded in the read name as: 
                                      @index1_index2_index3:original_read_name
                               output: demultiplexed fastq1 and fastq2 files 
          mapping: perform reads alignment
                             input: fastq files, separated by comma for each paired end
                             output: position sorted bam file, mapping qc stat and fragment.bed
          call_peak: call peaks for aggregated data
                               input: BAM file path
                               output: peaks in bed file format
          get_mtx: build raw peak by barcode matrix
                             input: features/peak file path
                             output: sparse matrix in Matrix Market format
          aggr_signal: generate aggregated signal, which can be upload to and view
                                 in genome browser
                                 input: require BAM file path
                                 output: bw and bedgraph file
          qc_per_barcode: quality control per barcode
                                    input: fragment.txt file and peak/feature file, separated by comma
                                    output: qc_per_barcode.summary
          process: processing data - including dex_fastq, mapping, call_peak, get_mtx,
                                aggr_signal, qc_per_barcode abd call_cell
                                input: fastq files for both reads and index, separated by comma like:
                                       fastq1,fastq2,index_fastq1,index_fastq2, index_fastq3...; 
                                output: cell peak matrix and all intermediate results 
          process_no_dex: processing data without demultiplexing
                                input: demultiplexed fastq files for both reads and index, separated by comma like:
                                       fastq1,fastq2; 
                                output: cell peak matrix and all intermediate results 
          process_with_bam: processing from bam file
                                input: bam file for aggregated data 
                                output: cell peak matrix and all intermediate results 
          call_cell: cell calling
                               input: raw peak barcode sparse matrix file path
                               output: filtered peak by cell matrix
          get_bam4Cells: extract bam file for cell barcodes and calculate mapping stats
                               input: bam file for all barcodes and a barcodes.txt file, separated by comma
                               output: bam file and mapping stats (optional) for cell barcodes                          
          clustering: cell clustering
                               input: filtered peak by cell matrix file path
                               output: seurat objects with clustering label in the metadata (.rds file) and 
                                       barcodes with cluster labels (bed file)
          motif_analysis: doing motif analysis
                               input: filtered peak by cell matrix file path
                               output: TF by cell matrix indicating TF accessibility (chromVAR object)
          runDA: doing differential accessibility analysis
                           input: seurat_obj.rds path from clustering analysis
                           output: differential peaks in txt format saved at the same directory as seurat_obj.rds
          runGO: doing GO analysis
                           input: result of runDA module (.txt file)
                           output: enriched GO terms in .xlsx saved at the same directory as the input file
          runCicero: run cicero for calculating gene activity score and predicting interactions
                           input: seurat_obj.rds path from clustering analysis
                           output: gene activity in .rds format and predicted interactions in .txt format

          split_bam: split bam file into different clusters
                               input: barcodes with cluster label (.txt file, outputed from clustering)
                               output: .bam (saved under downstream/CELL_CALLER/data_by_cluster), .bw, .bedgraph (save under output/signal/) file for each cluster
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
          convert10xbam: convert bam in 10x style to bam in scATAC-pro style 
                         input: bam file (position sorted) in 10x style
                         output: position sorted bam file in scATAC-pro style, mapping qc stat and fragment.bed
          mergePeaks: merge different peaks (called from differnt samples or conditions) if the distance is
                            less than a given #basepairs (200 if not specified) 
                         input: peak files and a distance paramter separated by comma, like:
                                peakFile1,peakFile2,peakFile3,200
                         output: merged peaks saved in file output/peaks/merged.bed
          reconstMtx: reconstruct peak-cell matrix given peak file, fragments.txt file, and barcodes.txt file 
                         input: differnt files separated by comma, like:
                                peakFilePath,fragmentFilePath,barcodesPath
                         output: a reconstructed peak-cell matrix saved under the same path as barcodes.txt file
          visualize: interactively visualize the data through web browser
                         input: VisCello_obj directory (created by clustering module)
                         output: web browser pop up for interactively visualization"

       -i|--input INPUT : input data, different types of input data are required for different steps;
       -c|--conf CONFIG : configuration file for parameters (if exists) for each step
       [-o|--output_dir : folder to save results, default output/ under the curret directory; sub-folder will be created automatically for each step
       [-h|--help]: help
       [-v|--version]: version




Run through docker or singularity
----------------------------------
In case you have problem in installing dependencies, you can run it without installing dependencies in **one of** following options:

1. Run the pre-built dockerized version [here](https://hub.docker.com/r/wbaopaul/scatac-pro) (which is not automatically updated)

2. Build your own docker image using the Dockfile in this repository (under scripts/install/Docker/, which is automatically updated) by running:

```
$ cd scripts/install/Docker
$ docker build -t YOU_PREFERED_NAME .
$ docker run -v YOUR_WORK_DIR:/software -it YOUR_PREFERED_NAME 
$ scATAC-pro --help


```
3. Run it through singularity (which is more fridenly with HPC and linux server) by running:

```
$ singularity pull -F docker://wbaopaul/scatac-pro 
## will output scatac-pro_latest.sif

$ singularity run -H YOUR_WORK_DIR --cleanenv scatac-pro_latest.sif
$ scATAC-pro --help

```

4. To use it on HPC cluster:

```
# write a script mapping.sh for mapping as an example:
#!/bin/bash
module load singularity
singularity pull -F docker://wbaopaul/scatac-pro
## will output scatac-pro_latest.sif
singularity exec -H YOUR_WORK_DIR --cleanenv scatac-pro_latest.sif scATAC-pro -s mapping
                            -i fastq_file1,fastq_file2,fastq_file3 -c configure_user.txt
# and then qsub mapping.sh


```

- **NOTE**: YOUR_WORK_DIR is your working directory, where the outputs will be saved and all data under YOUR_WORK_DIR will
be available to scATAC-pro

- **NOTE**: all inputs including data paths specified in configure_user.txt should be available under YOUR_WORK_DIR

FAQs
--------------
- [How to proceed using 10x cellranger-atac output?](https://github.com/wbaopaul/scATAC-pro/wiki/FAQs)
- [How to merge differnt peaks called from different sampels or conditions?](https://github.com/wbaopaul/scATAC-pro/wiki/FAQs)
- [How to reconstruct peak-by-cell matrix after updating peak file?](https://github.com/wbaopaul/scATAC-pro/wiki/FAQs)



Citation
--------------------------------------
Yu W, Uzun Y, Zhu Q, Chen C, Tan K. *scATAC-pro: a comprehensive workbench for single-cell chromatin accessibility sequencing data.* bioRxiv.org; 2019 
doi: https://doi.org/10.1101/824326 
