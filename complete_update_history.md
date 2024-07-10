## Complete Update History
-Version 1.5.2
    * Be compatible with Seurat v5
    * Add *report_dynamic* module: Interactive report with shiny features (still in testing)
    * *call_peak*: change to skip peak extension;
- Version 1.5.1 
    * *Integrate* module takes SampleSheet.csv file as input, where sample names, paths of peaks, fragments, and cell barcodes can be specified
    * Enable 4/5 bp shift for each read by setting SHIFT_READS_IN_BAM as TRUE in configure_user.txt file
    * Correctted an error for single-end sequencing data in the mapping module
- Version 1.5.0 released
    * Donot annotate peaks in the seurat object (peaks are still annotated for modules *runDA*, *runGO* and *visualize*)
    * Add new module *process_from_align* to do processing from alignment step (including aligment), given trimmed demultiplexed fastq files
    * Enable to change default expected doublet rate in the configure_user.txt file
- Version 1.4.4 released
    * Only consider standard chromosomes in the *qc_per_barcode* module 
    * Correct a minor bug in the *qc_per_barcode* module
    * Add version# in the html report
    * Clean and correct a minor bug in the *trimming* module
- Version 1.4.3 released
    * add new module *reprocess_cellranger_output* to reprocess scATAC-seq data originally processed by cellranger 
- Version 1.4.2 released
    * correct bugs in processing single-end sequencing data
    * correct minor underestimate of overall FrIP in the qc *report*
    * correct a bug in *process_with_bam* module
    * correct a bug in *integrate* module
- Version 1.4.1 released
    * *qc_per_barcode*: add tss enrichment score per cell into the QC metrics
    * update tutorial
    * correct a bug relate to input filepath for report module 
    * record input file paths for the *integration* module
    * using BAMPE for macs2 for paired-end sequencing
- Version 1.4.0 released 
    * new module added: labelTransfer (for cell annotation) from scRNA-seq
- Version 1.3.1 released
    * *rmDoublets*: a new module added, to remove doublets
    * *clustering*: accepts seurat obj (in .rds format) as input as well
- Version 1.3.0 released
    * *qc_per_barcode*: add tss enrichment score per cell into the QC metrics
    * *call_cell*: enable filtering barcodes by tss enrichment score
    * fragments file indexed by tabix (named fragments.tsv.gz)
    * *footprint* module: suppoort comparison of any two sets of cell clusters)
    * *motif_analysis* and *runDA*: accept seurat object in .rds format as input
    * *integrate*: rename cell name for each sample to avoid shared barcodes among samples; enable a distance parameter to merge peaks
    * *integrate_mtx*: added, as an alias of previous *integrate_seu* module
    * *report*: rearranged some plots and enabled output cicero interaction plot for a specific gene (specify it through **Cicero_Plot_Region* parameter in the configure_user.txt file)

- Version 1.2.1 released
    * new module *addCB2bam*: add cell barcode (CB) tag to a give bam file, new bam file will be saved in the same folder as the input bam (with name *_withCBtag.bam)
    * save .rds file for matrix and correct bug of calculating insert size 

- Version: 1.2.0 released:
    * update footprint dependency *rgt-hint* module to python3
    * save qc statistics in html report into tables, and peak calling summary inf added in the report
    * add qc per cell to seurat obj metadata as: total.unique.frags, frac.peak, frac.mito,
      frac.tss, frac.promoter, and frac.enhancer
- VERSION 1.1.4 released
    * *demplx_fastq*: the input supports directory path of 10x fastq files
- VERSION: 1.1.3 released
    * *runGO*: update background genes to be all genes associated with any peak
- May, 2020 --VERSION **1.1.2** released
    * *integrate*: add VFACS (Variable Features Across ClusterS) option for the integration module,
      **which reselect variable features across cell clusters after an initial clustering, followed by 
        another round of dimension reduction and clustering**, specify *Integrate_by = VFACS* in configure file
    * *clustering*: filter peaks before clustering (accessible in less than 0.5% of cells) and
       remove rare peaks (accessible in less than 1% of cells) from the variable features list
    * *reConsMtx*: enable specifying a path for saving reconstructed matrix (optional)
- VERSION **1.1.1** released
- March,April, 2020
    * *get_mtx*: it requires two input files: a fragments.txt file and a peak file, separated by comma
    * annotate peak as overlapped with a gene Tss if the corresponding distance <= 1000bp; mark peak with a gene if their distance <= 100kb
    * update DA, fix bug of using covariate
    * using mefa4::Melt instead of melt -- better for large sparse matrix
    * add PEAK_CALLER prefix to qc_per_barcode.txt filename
    * fix a bug of file location of tmpJob when calling cells
- VERSION **1.1.0** released
- Feb, 2020
    * *integrate* module enables 3 options: seurat, harmony and pool
    * new module *visualize*, allowing interactively explore and analyze the data
    * *footprint* module supports one-vs-rest comparison and provides result in heatmap
    * module *runDA* changed to use group name as the input (e.g. "0:1,2" or "one,rest")
    * installed rgt-hint (for footprinting analysis) using miniconda3
    * added module *process_with_bam*, allowing processing from aggregated bam file
    * enabled data integration from peaks files, assuming all data sets are processed using scATAC-pro. Output matrix with the same merged peaks/features and the previously called cells, along with an integrated seurat object
    * added new parameters in the configuration file: Top_Variable_Features, REDUCTION, nREDUCTION
    * enabled all clustering methods mentioned in the manuscript, along with kmeans clustering on principal components
    * file path changed to like downstreame_analysis/PEAK_CALLER/CELL_CALLER/..., indicating peak caller
- Jan, 2020
    * added a new module *mergePeaks* to merge different peak files called from different data sets
    * added a new module *reConstMtx* to reconstruct peak-by-cell matrix given a peak file, a fragment file and a barcodes.txt file
- Dec, 2019
    * corrected an error due to using older version of chromVAR
    * corrected a bug for demultiplexing multiple index files
    * added a module *convert10xbam* to convert 10x position sorted bam file to scATAC-pro file format
    * updated module *get_bam4Cells*, with the inputs as a bam file and a txt file of barcodes, separated by comma

