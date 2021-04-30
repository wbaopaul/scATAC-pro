source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local('dsAnalysis_utilities.R')

args = commandArgs(T)
seuratObj_file = args[1]
drate = as.numeric(args[2])
if(is.null(drate)) drate = 0.4
seurat.obj = readRDS(seuratObj_file)

seurat.obj = FindDoublets_Atac(seurat.obj, PCs = 1:10, exp_rate = drate,
             sct = F)

## plot a umap with Doublet/Singlet
p0 <- DimPlot(seurat.obj, group.by = 'Doublet_Singlet')
output_dir = dirname(seuratObj_file)
output_seu_file = gsub('.rds', '_doubletsRemoved.rds', seuratObj_file,
 fixed = T)
ggsave(p0, file = paste0(output_dir, '/umap_singlets_doublets.eps'),
      device = 'eps', width = 6, height = 6)

seurat.obj <- subset(seurat.obj, Doublet_Singlet == 'Singlet')
saveRDS(seurat.obj, file = output_seu_file)
