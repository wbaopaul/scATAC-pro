source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local('dsAnalysis_utilities.R')

library(edgeR)
args = commandArgs(T)
seuratObj_file = args[1]
output_dir = args[2]
group1 = args[3]
group2 = args[4]
test_use = args[5]

seurat.obj = readRDS(seuratObj_file)

seurat.obj$active_clusters = as.character(seurat.obj$active_clusters)
Idents(seurat.obj) <- seurat.obj$active_clusters

group1 = tolower(group1)
group2 = tolower(group2)

confVar = NULL
slot = 'data'
if(test_use %in% c('DESeq2', 'LR', 'negbinom')){
  slot = 'counts'
  confVar = 'nCount_ATAC'
}

cls = sort(unique(seurat.obj$active_clusters))
group1 = unlist(strsplit(group1, ':'))
group2 = unlist(strsplit(group2, ':'))
cn = colnames(seurat.obj)

if(group1[1] == 'one') {
   markers = NULL
   for(cluster0 in cls){
        cells1 = cn[which(seurat.obj$active_clusters == cluster0)]
        if(group2[1] == 'rest') {
          cells2 = cn[which(seurat.obj$active_clusters != cluster0)]
          id2 = NULL
        }else{
          cells2 = cn[which(seurat.obj$active_clusters %in% group2)]
          id2 = group2
        }
        if(length(cells1) <= 10 || length(cells2) <= 10) next
        mm = FindMarkers(seurat.obj, slot = slot, 
                         ident.1 = cluster0, ident.2 = id2,
                         test.use = test_use, logfc.threshold = 0.0, 
                         max.cells.per.ident = 500, 
                         only.pos = T, latent.vars = confVar)

        mm$cluster = cluster0
        mm$fdr = p.adjust(mm$p_val, method = 'fdr')
        mm$peak = rownames(mm)
        markers = rbind(markers, mm)

   }
}else{
  cells1 = cn[which(seurat.obj$active_clusters %in% group1)]
  if(length(cells1) <= 10) stop('Not enough cells in group1')
    if(group2[1] == 'rest'){
        cells2 = cn[which(!seurat.obj$active_clusters %in% group1)]
        id2 = setdiff(cls, group1)
        markers = FindMarkers(seurat.obj, slot = slot,
                              ident.1 = group1, ident.2 = id2, test.use = test_use, 
                              logfc.threshold = 0.0, max.cells.per.ident = 500,
                              only.pos = F, latent.vars = confVar)
        markers$cluster = ifelse(markers$avg_logFC > 0, group1, group2)
        markers$fdr = p.adjust(markers$p_val, method = 'fdr')
    }else{
        cells2 = cn[which(seurat.obj$active_clusters %in% group2)]
        if(length(cells2) <= 10) stop('Not enough cells in group2')
        markers = FindMarkers(seurat.obj, slot = slot, 
                              ident.1 = group1, ident.2 = group2, test.use = test_use,
                              max.cells.per.ident = 500, logfc.threshold = 0.0, 
                              only.pos = F, latent.vars = confVar)
        markers$cluster = ifelse(markers$avg_logFC > 0, group1, group2)
        markers$fdr = p.adjust(markers$p_val, method = 'fdr')
    }
  
  markers$peak = rownames(markers)
}


markers = data.table(markers)
markers[, 'peak0' := unlist(strsplit(peak, ','))[1], by = peak]
markers[, 'chr' := unlist(strsplit(peak0, '-'))[1], by = peak0]
markers[, 'start' := unlist(strsplit(peak0, '-'))[2], by = peak0]
markers[, 'end' := unlist(strsplit(peak0, '-'))[3], by = peak0]

setcolorder(markers, c('chr', 'start', 'end', 'p_val','avg_logFC','pct.1','pct.2', 
                       'p_val_adj', 'fdr', 'cluster', 'peak', 'peak0'))

markers = markers[fdr <= 0.05, ]
write.table(markers, file = paste0(output_dir, '/differential_accessible_features_', args[3], '_vs_', args[4], '.txt'), sep = '\t',
            quote = F, row.names = F)

