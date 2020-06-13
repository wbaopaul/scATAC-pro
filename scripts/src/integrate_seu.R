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

message(paste('Integrate by', integrate_by))

mtx_files = unlist(strsplit(mtx_files, ','))

tss_ann <- fread(tss_path, header = F)
#names(tss_ann)[c(1:4,7)] <- c('chr', 'start', 'end', 'gene_name', 'gene_type')
#tss_ann <- tss_ann[gene_type %in% c('miRNA', 'lincRNA', 'protein_coding'), ]
names(tss_ann)[c(1:4)] <- c('chr', 'start', 'end', 'gene_name')

len = length(mtx_files)
mtx_list = list()
for(i in 1:len){
    file0 = mtx_files[i]
    mtx = read_mtx_scATACpro(file0)
    mtx = assignGene2Peak(mtx, tss_ann)
    mtx_list[[i]] = mtx
}
seurat.obj <- run_integration(mtx_list, integrate_by = integreate_by,
                            top_variable_features = top_variable_features, 
                            norm_by = norm_by, nREDUCTION = nREDUCTION,
                            minFrac_in_cell = 0.01, min_depth = 1000,
                            max_depth = 50000, reg.var = 'nCount_ATAC',
                            anchor.features = NULL,
                            resolution = 0.6)


saveRDS(seurat.obj, file = paste0(output_dir, '/seurat_obj_', integrate_by, '.rds'))

#output bacrcode cluster information
bc_cls = data.table('Barcode' = rownames(seurat.obj@meta.data), 'Cluster' = seurat.obj@meta.data$active_clusters)
setkey(bc_cls, Cluster)

write.table(bc_cls, file = paste0(output_dir, '/cell_cluster_table_', integrate_by, '.txt'), sep = '\t',
            quote = F, row.names = F)

cg <- DimPlot(seurat.obj, reduction = 'umap', group.by = 'active_clusters', label = T) + 
  theme(legend.text = element_text(size = 17)) 
if(length(unique(seurat.obj$active_clusters)) < 10) cg = cg + scale_color_brewer(palette = "Paired")

cg1 <- DimPlot(seurat.obj, reduction = 'umap', group.by = 'sample') + 
  theme(legend.text = element_text(size = 17)) + scale_color_brewer(palette = "Set1")

pfname = paste0(output_dir, '/umap_clusters_', integrate_by, '.eps')
   ggsave(CombinePlots(plots = list(cg1, cg)), file = pfname, device = 'eps', width = 14, height = 6)

