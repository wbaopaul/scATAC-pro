source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local('dsAnalysis_utilities.R')

args = commandArgs(T)
inputSeurat_atac = args[1] ## a seurat.rds file, with pca conducted and metadata Cell_Type 
inputSeurat_rna = args[2] ## a seurat.rds file, with pca conducted
GENOME_NAME = args[3]
gene_gtf_file = args[4]

## if gtf file not provided, using R bioconductor packages for gene annotation
if(!file.exists(gene_gtf_file)){
  if(!grepl(GENOME_NAME, pattern = 'hg19|hg38|mm10|mm9', ignore.case = T)){
    stop('Genome does not belong to any of hg19,hg38,mm9 or mm10, 
         please provide .gtf file for gene annotation!')
  }
  if(grepl(GENOME_NAME, pattern = 'mm10', ignore.case = T)) {
    library(EnsDb.Mmusculus.v79)
    ens.ann <- EnsDb.Mmusculus.v79
  }
  if(grepl(GENOME_NAME, pattern = 'mm9', ignore.case = T)) {
    library(EnsDb.Mmusculus.v75)
    ens.ann <- EnsDb.Mmusculus.v75
  }
  
  if(grepl(GENOME_NAME, pattern = 'hg38', ignore.case = T)) {
    library(EnsDb.Hsapiens.v86)
    ens.ann <- EnsDb.Hsapiens.v86
  }
  
  if(grepl(GENOME_NAME, pattern = 'hg19', ignore.case = T)) {
    library(EnsDb.Hsapiens.v75)
    ens.ann <- EnsDb.Hsapiens.v75
  }
  
  gene_ann <- ensembldb::genes(ens.ann)
  gene_ann <- keepStandardChromosomes(gene_ann, pruning.mode = 'coarse')
  gene_ann = data.frame(gene_ann)
  gene_ann = subset(gene_ann, gene_biotype %in% c('protein_coding', 'miRNA', 'lincRNA'),
                    select = c('seqnames', 'start', 'end',
                               'strand', 'gene_biotype', 'gene_name'))
  
  names(gene_ann)[1:3] = c('chr', 'gene_start', 'gene_end')
  gene_ann = data.table(gene_ann)
  gene_ann = gene_ann[!duplicated(gene_name)]
  gene_ann$chr = paste0('chr', gene_ann$chr)
  gene_ann = subset(gene_ann, select = c('chr', 'gene_start', 'gene_end',
                                         'strand', 'gene_name'))
  
  gene_ann[chr == 'chrMT']$chr = 'chrM'
  
}else{
  gene_ann = fread(gene_gtf_file, sep = '\t')
  gene_ann = gene_ann[V3 == 'gene']
  gene_ann[, 'gene_name' := unlist(strsplit(V9, ';'))[3], by = V9]
  gene_ann[, 'gene_name' := gsub("\"", "", gene_name), by = gene_name]
  gene_ann[, 'gene_name' := unlist(strsplit(gene_name, ' '))[3], by = gene_name]
  names(gene_ann)[1] = 'chr'
  gene_ann = subset(gene_ann, select = c(chr, V4, V5, V7, gene_name))
  names(gene_ann)[2:4] = c('start', 'end', 'strand')
  chrs = standardChromosomes(makeGRangesFromDataFrame(gene_ann))
  gene_ann = gene_ann[chr %in% chrs]
  gene_ann = gene_ann[!duplicated(gene_name)]
  names(gene_ann)[2:3] = c('gene_start', 'gene_end')
}

seurat.rna = readRDS(inputSeurat_rna)
seurat.atac = readRDS(inputSeurat_atac)

seurat.atac = labelTransfer_R(seurat.atac, seurat.rna, gene_ann,
                              rna_ann_var = 'Cell_Type', include_genebody = T)

p1 <- DimPlot(seurat.atac, group.by = "Predicted_Cell_Type",
              label = TRUE, repel = TRUE) + ggtitle("scATAC-seq") 

p2 <- DimPlot(seurat.rna, group.by = "Cell_Type",
              label = TRUE, repel = TRUE) + ggtitle("scRNA-seq") 

outputPath = dirname(inputSeurat_atac)
outputName = basename(inputSeurat_atac)
saveRDS(seurat.atac, file = paste0(outputPath, '/updated_', outputName))

pcomb = gridExtra::grid.arrange(p1, p2, nrow = 1)

ggsave(pcomb, filename = paste0(outputPath, '/umap_with_predicted_cell_type.eps'),
       device = 'eps', height = 6, width = 13)


