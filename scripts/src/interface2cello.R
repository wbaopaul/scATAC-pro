source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local('dsAnalysis_utilities.R')

library(ggplot2)
library(parallel)
args = commandArgs(T)

seuratPath = args[1]
assay4cello = args[2]
output_dir = dirname(seuratPath)

seurat.atac = readRDS(seuratPath)
rr = rownames(seurat.atac)
sele.features = VariableFeatures(seurat.atac)
rr1 = rr[grepl(rr, pattern = 'Tss')]
sele.features = unique(c(sele.features, rr1))

inputs = prepInput4Cello(seurat.atac@assays$ATAC@counts[rr %in% sele.features, ], 
                         seurat.obj = seurat.atac,
                         norm_mtx = seurat.atac@assays$ATAC@data[rr %in% sele.features, ],
                         cello.name = 'scATAC_withGene2Peak',
                         assay = assay4cello, extraDims = c(100, 80, 50, 30, 20),
                         vars.to.regOnPca = 'nCount_ATAC',
                         vFeatures = NULL)


celloPath = paste0(output_dir, '/VisCello_obj')
system(paste0('mkdir -p ', celloPath))
saveRDS(inputs$eset, file = paste0(celloPath, '/eset.rds'))
saveRDS(inputs$clist, file = paste0(celloPath, '/clist.rds'))



