source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local('dsAnalysis_utilities.R')

args = commandArgs(T)
inputObj_file = args[1] ## a mtx.rds or a seurat.rds file
drate = as.numeric(args[2])
if(is.null(drate)) drate = 0.03
input.obj = readRDS(inputObj_file)
if(any(class(input.obj) == 'Seurat')) {
    seurat.obj = input.obj
}else{
    #the input is mtx, create a seurat obj then
   seurat.obj = runSeurat_Atac(input.obj, npc = 30, norm_by = 'tf-idf',
                               top_variable_features = 5000, 
                               reg.var = 'nCount_ATAC') 
   seurat.obj = FindNeighbors(seurat.obj, reduction = 'pca', dims = 1:30, k.param = 20)    
   seurat.obj = FindClusters(seurat.obj, resolusion = 0.6)    
   seurat.obj = RunTSNE(seurat.obj, dims = 1:30, reduction = 'pca', check_duplicates = FALSE)
   seurat.obj = RunUMAP(seurat.obj, dims = 1:30, reduction = 'pca', verbose = F)

}

seurat.obj = FindDoublets_Atac(seurat.obj, PCs = 1:10, exp_rate = drate,
             sct = F)

## plot a umap with Doublet/Singlet
p0 <- DimPlot(seurat.obj, group.by = 'Doublet_Singlet')
output_dir = dirname(inputObj_file)
ggsave(p0, file = paste0(output_dir, '/umap_with_doublets.eps'),
      device = 'eps', width = 6, height = 6)
saveRDS(seurat.obj, file = paste0(output_dir, '/seurat_obj_with_doublets.rds'))

seurat.obj <- subset(seurat.obj, Doublet_Singlet == 'Singlet')

output_seu_file = paste0(output_dir, '/seurat_obj_doubletsRemoved.rds') 
output_mtx_file = paste0(output_dir, '/matrix_doubletsRemoved.rds') 
output_barcode_file = paste0(output_dir, '/barcodes_doubletsRemoved.txt') 
saveRDS(seurat.obj@assays$ATAC@counts, file = output_mtx_file)
write.table(colnames(seurat.obj), file = output_barcode_file,
            quote = F, row.names = F, col.names = F, sep = '\t')
saveRDS(seurat.obj, file = output_seu_file)

