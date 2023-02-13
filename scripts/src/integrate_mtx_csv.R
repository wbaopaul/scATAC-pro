source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local('dsAnalysis_utilities.R')

library(ggplot2)
library(parallel)
args = commandArgs(T)
mtx_files = args[1]
k = (args[2])
output_dir = args[3]
genome_name = args[4]
tss_path = args[5]
norm_by = args[6]
REDUCTION = args[7]
nREDUCTION = as.integer(args[8])
top_variable_features = as.numeric(args[9])
integrate_by = args[10]
sampleSheet = args[11]

message(paste('Integrate by', integrate_by))

mtx_files = unlist(strsplit(mtx_files, ','))

tss_ann <- fread(tss_path, header = F)
names(tss_ann)[c(1:4)] <- c('chr', 'start', 'end', 'gene_name')

sample_inf = fread(sampleSheet)
sampleNames = sample_inf$sample_name

len = length(mtx_files)
seurat_list = list()
for(i in 1:len){
    file0 = mtx_files[i]
    sample0 = sampleNames[i]
    mtx = readRDS(file0)
    message(paste0('Working on processing ', sample0, ' ...'))
    ## rename cell names in case of shared barcodes among samples
    colnames(mtx) = paste0(sample0, '_', colnames(mtx)) 
    seurat0 = runSeurat_Atac(mtx, npc = nREDUCTION, norm_by = norm_by, 
                                top_variable_features = top_variable_features, 
                                reg.var = 'nFeature_ATAC')
    seurat0$sampleName = sample0
    seurat_list[[sample0]] = seurat0
}

message(paste0('Working on integration ...'))
seurat.obj <- run_integrateSeuObj(seurat_list, integrate_by = integrate_by,
                            top_variable_features = top_variable_features, 
                            norm_by = norm_by, nREDUCTION = nREDUCTION,
                            reg.var = 'nFeature_ATAC',
                            resolution = 0.6)

saveRDS(seurat.obj, file = paste0(output_dir, '/seurat_obj_', integrate_by, '.rds'))

#output bacrcode cluster information
bc_cls = data.table('Barcode' = rownames(seurat.obj@meta.data), 'Cluster' = seurat.obj@meta.data$active_clusters)
setkey(bc_cls, Cluster)

write.table(bc_cls, file = paste0(output_dir, '/cell_cluster_table_', integrate_by, '.tsv'), sep = '\t',
            quote = F, row.names = F)

getPalette1 = colorRampPalette(brewer.pal(9, "Paired"))
myColors1 = getPalette1(len)
getPalette2 = colorRampPalette(brewer.pal(9, "Set1"))
seurat.obj$active_clusters = as.character(seurat.obj$active_clusters)
myColors2 = getPalette2(length(unique(seurat.obj$active_clusters)))

cg1 <- DimPlot(seurat.obj, reduction = 'umap', group.by = 'sampleName') + 
  theme(legend.text = element_text(size = 17)) + scale_color_manual(values = myColors1)
cg2 <- DimPlot(seurat.obj, reduction = 'umap', group.by = 'active_clusters', label = T) + 
  theme(legend.text = element_text(size = 17)) + scale_color_manual(values = myColors2) 


pcomb = gridExtra::grid.arrange(cg1, cg2, nrow = 1)
pfname = paste0(output_dir, '/umap_clusters_', integrate_by, '.eps')

#ggsave(CombinePlots(plots = list(cg1, cg)), file = pfname, device = 'eps', width = 14, height = 6)
ggsave(pcomb, file = pfname, device = 'eps', width = 14, height = 6)

