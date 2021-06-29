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
    stop('Genome is not belong to any of hg19,hg38,mm9 or mm10, 
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
  chrs = 1:22
  chrs = c(chrs, 'X', 'Y', 'M')
  gene_ann = gene_ann[chr %in% chrs]
  gene_ann = gene_ann[!duplicated(gene_name)]
  names(gene_ann)[2:4] = c('gene_start', 'gene_end', 'strand')
  gene_ann[, 'chr' := paste0('chr', chr)]
  
}


if(F){
  ## download gtf file if not provided
  if(!file.exists(gene_gtf_file)){
    print('gene annotation gtf file not provided, I will try to download one:')
    err = 0
    if(grepl(GENOME_NAME, pattern = 'mm9', ignore.case = T)) {
      err <- tryCatch(download.file('ftp://ftp.ensembl.org/pub/release-67/gtf/mus_musculus/Mus_musculus.NCBIM37.67.gtf.gz', temp),
                      error = function(e) {
                        print("Cannot download ftp://ftp.ensembl.org/pub/release-67/gtf/mus_musculus/Mus_musculus.NCBIM37.67.gtf.gz!")
                        return(1)})
    }
    if(grepl(GENOME_NAME, pattern = 'mm10', ignore.case = T)) {
      err <- tryCatch(download.file('ftp://ftp.ensembl.org/pub/release-95/gtf/mus_musculus/Mus_musculus.GRCm38.95.gtf.gz', temp),
                             error = function(e) {
                               print("Cannot download ftp://ftp.ensembl.org/pub/release-95/gtf/mus_musculus/Mus_musculus.GRCm38.95.gtf.gz!")
                               return(1)})
    }
    if(grepl(GENOME_NAME, pattern = 'hg38', ignore.case = T)) {
      err <- tryCatch(download.file('ftp://ftp.ensembl.org/pub/release-95/gtf/homo_sapiens/Homo_sapiens.GRCh38.95.gtf.gz', temp),
                      error = function(e) {
                        print("Cannot download ftp://ftp.ensembl.org/pub/release-95/gtf/homo_sapiens/Homo_sapiens.GRCh38.95.gtf.gz!")
                        return(1)})
    }
    
    if(grepl(GENOME_NAME, pattern = 'hg19', ignore.case = T)) {
      err <- tryCatch(download.file('ftp://ftp.ensembl.org/pub/release-67/gtf/homo_sapiens/Homo_sapiens.GRCh37.67.gtf.gz', temp),
                      error = function(e) {
                        print("Cannot download ftp://ftp.ensembl.org/pub/release-67/gtf/homo_sapiens/Homo_sapiens.GRCh37.67.gtf.gz!")
                        return(1)})
    }
    
    if(err == 1) stop('Download failed! Please provide a gtf file to run this module!')
  }
  
  gene_ann = fread(gene_gtf_file, sep = '\t')
  gene_ann = gene_ann[V3 == 'gene']
  gene_ann[, 'gene_name' := unlist(strsplit(V9, ';'))[3], by = V9]
  gene_ann[, 'gene_name' := gsub("\"", "", gene_name), by = gene_name]
  gene_ann[, 'gene_name' := unlist(strsplit(gene_name, ' '))[3], by = gene_name]
  names(gene_ann)[1] = 'chr'
  gene_ann = subset(gene_ann, select = c(chr, V4, V5, V7, gene_name))
  chrs = 1:22
  chrs = c(chrs, 'X', 'Y')
  gene_ann = gene_ann[chr %in% chrs]
  gene_ann = gene_ann[!duplicated(gene_name)]
  names(gene_ann)[2:4] = c('gene_start', 'gene_end', 'strand')
  gene_ann[, 'chr' := paste0('chr', chr)]
  
}

seurat.rna = readRDS(inputSeurat_rna)
seurat.atac = readRDS(inputSeurat_atac)

atac.mtx = seurat.atac@assays$ATAC@counts
rn = rownames(atac.mtx)
rownames(atac.mtx) <- sapply(rn, function(x) unlist(strsplit(x, ','))[1])
activity.matrix = generate_gene_cisActivity(gene_ann = gene_ann,
                                            atac.mtx, 
                                            include_body = T)

seurat.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)

DefaultAssay(seurat.atac) <- "ACTIVITY"
seurat.atac <- NormalizeData(seurat.atac)
seurat.atac <- ScaleData(seurat.atac, features = rownames(seurat.atac))
seurat.atac <- FindVariableFeatures(seurat.atac)

DefaultAssay(seurat.atac) <- "ATAC"
seurat.atac$tech = 'ATAC'
seurat.rna$tech = 'RNA'

## transfer label 
genes4anchors = VariableFeatures(object = seurat.rna)
#genes4anchors = NULL
transfer.anchors <- FindTransferAnchors(reference = seurat.rna,
                                        query = seurat.atac,
                                        features = genes4anchors,
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 5)


celltype.predictions <- TransferData(anchorset = transfer.anchors,
                                     refdata = seurat.rna$Cell_Type,
                                     weight.reduction = seurat.atac[["pca"]],
                                     dims = 1:ncol(seurat.atac[["pca"]]),
                                     k.weight = 50)
celltype.predictions = subset(celltype.predictions, select = c('predicted.id', 'prediction.score.max'))
names(celltype.predictions)[1] = 'Predicted_Cell_Type'
seurat.atac <- AddMetaData(seurat.atac, metadata = celltype.predictions)
rm(transfer.anchors)

p1 <- DimPlot(seurat.atac, group.by = "Predicted_Cell_Type",
              label = TRUE, repel = TRUE) + ggtitle("scATAC-seq") 

p2 <- DimPlot(seurat.rna, group.by = "Cell_Type",
              label = TRUE, repel = TRUE) + ggtitle("scRNA-seq") 

seurat.atac[["ACTIVITY"]] <- NULL ## don't save activity assay

outputPath = dirname(inputSeurat_atac)
outputName = basename(inputSeurat_atac)
saveRDS(seurat.atac, file = paste0(outputPath, '/updated_', outputName))

pcomb = gridExtra::grid.arrange(p1, p2, nrow = 1)

ggsave(pcomb, filename = paste0(outputPath, '/umap_with_predicted_cell_type.eps'),
       device = 'eps', height = 6, width = 13)


