source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local('dsAnalysis_utilities.R')

library(parallel)
args = commandArgs(T)
cluster1 = args[1]
cluster2 = args[2]
test_use = args[3]
output_dir = args[4]


if(file.exists(paste0(output_dir, '/seurat_obj_withCluster.rds'))){
  seurat.obj = readRDS(paste0(output_dir, '/seurat_obj_withCluster.rds'))
}else{
  stop('Should do clustering analysis first!')
}

if(cluster2 == 'others') cluster2 = NULL
if(cluster1 == 'all') {
   cls = unique(seurat.obj$active_clusters)
   markers = NULL
   if(cluster2 == 'others' || cluster2 == 'all') cluster2 = NULL
   for(cluster0 in cls){
        mm = FindMarkers(seurat.obj, group.by = active_clusters, ident.1 = cluster0, ident.2 = cluster2,
                test.use = test_use, max.cell.per.ident = 500, only.pos = T)
        markers = rbind(markers, mm)
   }
}else{
    markers = FindMarkers(seurat.obj, group.by = active_clusters, ident.1 = cluster1, ident.2 = cluster2, 
    test.use = test_use, max.cell.per.ident = 500, only.pos = T)
}

write.table(markers, file = paste0(output_dir, '/differential_accessible_peaks.txt'), sep = '\t',
            quote = F, row.names = F)

