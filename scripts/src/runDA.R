source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local('dsAnalysis_utilities.R')


args = commandArgs(T)
seuratObj_file = args[1]
output_dir = args[2]
group1 = args[3]
group2 = args[4]
test_use = args[5]

seurat.obj = readRDS(seuratObj_file)

confVar = 'nCount_ATAC'
if(test_use == 'wilxon' || test_use == 'DESeq2') confVar = NULL

if(group2 == 'rest') group2 = NULL
if(group1 == 'all') {
   cls = unique(seurat.obj$active_clusters)
   markers = NULL
   for(cluster0 in cls){
        mm = FindMarkers(seurat.obj, group.by = "active_clusters", ident.1 = cluster0, ident.2 = group2,
                test.use = test_use, logfc.threshold = 0.0, max.cell.per.ident = 500, only.pos = T, latent.vars = confVar)

        mm$cluster = cluster0
        
        mm$peak = rownames(mm)
        markers = rbind(markers, mm)

   }
}else{
    if(is.null(group2)){
        markers = FindMarkers(seurat.obj, group.by = "active_clusters", ident.1 = group1, ident.2 = group2, 
                          test.use = test_use, logfc.threshold = 0.0, max.cell.per.ident = 500, only.pos = T, latent.vars = confVar)
        markers$cluster = group1
    }else{
        markers = FindMarkers(seurat.obj, group.by = "active_clusters", ident.1 = group1, ident.2 = group2, test.use = test_use, max.cell.per.ident = 500, logfc.threshold = 0.0, only.pos = F, latent.vars = confVar)
        markers$cluster = ifelse(markers$avg_logFC > 0, group1, group2)
    }
  
  markers$peak = rownames(markers)
}


markers = data.table(markers)
markers[, 'peak0' := unlist(strsplit(peak, ','))[1], by = peak]
markers[, 'chr' := unlist(strsplit(peak, '-'))[1], by = peak0]
markers[, 'start' := unlist(strsplit(peak, '-'))[2], by = peak0]
markers[, 'end' := unlist(strsplit(peak, '-'))[3], by = peak0]

setcolorder(markers, c('chr', 'start', 'end', 'p_val','avg_logFC','pct.1','pct.2', 'p_val_adj','cluster', 'peak'))
markers = markers[abs(avg_logFC) > 0.1, ]
write.table(markers, file = paste0(output_dir, '/differential_peak_cluster_table.txt'), sep = '\t',
            quote = F, row.names = F)

