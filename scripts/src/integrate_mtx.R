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
    if(grepl(file0, pattern = '.rds', fix = T)){
        mtx = readRDS(file0)
    }else{
        mtx = read_mtx_scATACpro(file0)
    }
#   mtx = assignGene2Peak(mtx, tss_ann)

    ## rename cell names in case of shared barcodes among samples
    colnames(mtx) = paste0('sample', i, '_', colnames(mtx))    

    mtx_list[[i]] = mtx
}
names(mtx_list) = paste0('sample', 1:len)
seurat.obj <- run_integration(mtx_list, integrate_by = integrate_by,
                            top_variable_features = top_variable_features, 
                            norm_by = norm_by, nREDUCTION = nREDUCTION,
                            minFrac_in_cell = 0.01, min_depth = 500,
                            max_depth = 50000, reg.var = 'nCount_ATAC',
                            resolution = 0.6)

# record input sample path as metadata to the seurat object
dpath = data.table('sample' = paste0('sample', 1:len),
                   'sample_path' = mtx_files)
setkey(dpath, sample)
seurat.obj$sample_path = dpath[J(seurat.obj$sample)]$sample_path

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

cg1 <- DimPlot(seurat.obj, reduction = 'umap', group.by = 'sample') + 
  theme(legend.text = element_text(size = 17)) + scale_color_manual(values = myColors1)
cg2 <- DimPlot(seurat.obj, reduction = 'umap', group.by = 'active_clusters', label = T) + 
  theme(legend.text = element_text(size = 17)) + scale_color_manual(values = myColors2) 


pcomb = gridExtra::grid.arrange(cg1, cg2, nrow = 1)
pfname = paste0(output_dir, '/umap_clusters_', integrate_by, '.eps')

#ggsave(CombinePlots(plots = list(cg1, cg)), file = pfname, device = 'eps', width = 14, height = 6)
ggsave(pcomb, file = pfname, device = 'eps', width = 14, height = 6)

