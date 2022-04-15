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
tss_path = args[3]  ## changed since 1.5.0

output_dir = dirname(seuratPath)

seurat.atac = readRDS(seuratPath)
rr = rownames(seurat.atac)
sele.features = VariableFeatures(seurat.atac)
mtx <- seurat.atac@assays$ATAC@counts[rr %in% sele.features, ]
mtx.norm <- seurat.atac@assays$ATAC@data[rr %in% sele.features, ]

## annote peaks with genes
tss_ann <- fread(tss_path, header = F)
names(tss_ann)[c(1:4)] <- c('chr', 'start', 'end', 'gene_name')
mtx = assignGene2Peak(mtx, tss_ann)
mtx.norm = assignGene2Peak(mtx.norm, tss_ann)
rr = rownames(mtx)
rr1 = rr[grepl(rr, pattern = 'Tss')]

inputs = prepInput4Cello(mtx[rr %in% rr1, ],
                         seurat.obj = seurat.atac,
                         norm_mtx = mtx.norm[rr %in% rr1, ], 
                         cello.name = 'scATAC_withGene2Peak',
                         assay = assay4cello, extraDims = c(100, 80, 50, 30, 20),
                         vars.to.regOnPca = 'nCount_ATAC',
                         vFeatures = NULL)


celloPath = paste0(output_dir, '/VisCello_obj')
system(paste0('mkdir -p ', celloPath))
saveRDS(inputs$eset, file = paste0(celloPath, '/eset.rds'))
saveRDS(inputs$clist, file = paste0(celloPath, '/clist.rds'))



