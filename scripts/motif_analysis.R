source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local('dsAnalysis_utilities.R')

library(parallel)
library(pheatmap)
library(viridis)

args = commandArgs(T)
mtx_file = args[1]
genome_name = args[2]
output_dir = args[3]

mtx = read_mtx_scATACpro(mtx_file)
mtx = filterMat(mtx)
genomeName = 'BSgenome.Hsapiens.UCSC.hg38'
if(grepl(genome_name, pattern = '38'))genomeName = 'BSgenome.Hsapiens.UCSC.hg38'
if(grepl(genome_name, pattern = '19'))genomeName = 'BSgenome.Hsapiens.UCSC.hg19'
if(grepl(genome_name, pattern = 'mm9'))genomeName = 'BSgenome.Mmusculus.UCSC.mm9'
if(grepl(genome_name, pattern = 'mm10'))genomeName = 'BSgenome.Mmusculus.UCSC.mm10'

ncore = detectCores()
obj = run_chromVAR(mtx, genomeName, max(1, ncore - 1))
saveRDS(obj, file = paste0(output_dir, '/chromVar_obj.rds'))

## load seurat object with cluster information
#seurat.obj = readRDS(paste0(output_dir, '/seurat_obj_withCluster.rds'))

## check enriched TFs for each cluster
variability <- computeVariability(obj)
variability = data.table(variability, stringsAsFactors = F)
variability = variability[order(-variability), ]

## plot enriched TFs in heatmap
sele.tfs = as.character(variability$name[1:50])
zscores = obj@assays$data$z
rnames = rownames(zscores)
rownames(zscores) = sapply(rnames, function(x) unlist(strsplit(x, '_'))[2])
sele.zscores = zscores[rownames(zscores) %in% sele.tfs, ]

metaData = seurat.obj@meta.data
metaData = data.table(metaData, keep.rownames = T)
setkey(metaData, seurat_cluster)

sele.zscores = sele.zscores[, metaData$rn]

ann_col = data.table(cluster = metaData$seurat_cluster)
colnames(ann_col) = 'cluster'
rownames(ann_col) = metaData$rn

ann_colors = list(metaData$seurat_cluster)
sele.zscores[sele.zscores > 10] = 10
sele.zscores[sele.zscores < -10] = -10

pdf(paste0(output_dir, '/top50_TFs_chromVAR_zscore.heatmap.pdf'),
    width = 8, height = 8)
pheatmap::pheatmap(sele.zscores, cluster_cols = F, show_colnames = F,
                   annotation_col = ann_col, color = viridis(100))
dev.off()



