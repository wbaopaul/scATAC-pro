source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local('dsAnalysis_utilities.R')

library(parallel)
args = commandArgs(T)
mtx_file = args[1]
cluster_method = args[2]
k = as.integer(args[3])
output_dir = args[4]

mtx = read_mtx_scATACpro(mtx_file)



## doing dimension reduction if it hasn't been done
if(!file.exists(paste0(output_dir, '/seurat_obj_withDR.rds'))){
  seurat.obj = doBasicSeurat(mtx, npc = 100, top.variable = 0.2)
  
  seurat.obj = RunTSNE(seurat.obj, dims = 1:50)
  seurat.obj = RunUMAP(seurat.obj, dims = 1:50, verbose = F)
  saveRDS(seurat.obj, file = paste0(output_dir, '/seurat_obj_withDR.rds'))
  
}else{
  seurat.obj = readRDS(paste0(output_dir, '/seurat_obj_withDR.rds'))
}

## clustering
if(cluster_method == 'seurat'){
  ## seurat implemented louvain algorithm
  seurat.obj = FindNeighbors(seurat.obj, reduction = 'pca')
  seurat.obj = FindClusters(seurat.obj, resolution = 0.2)
  seurat.obj$seurat_cluster = seurat.obj@active.ident
}

if(cluster_method == 'cisTopic'){
  nc = detectCores()
  cis.obj = run_cisTopic(mtx, nCores = nc - 1)
  sele.cisTopic <- selectModel(cis.obj, 
                               keepBinaryMatrix = F, keepModels = F)
  cell_topic <- t(modelMatSelection(sele.cisTopic, 'cell', 'Probability'))
  seurat.obj$cisTopic_cluster = generalCluster(cell_topic, method = 'hclust', k = k)
}

if(cluster_method == 'LSI'){
  seurat.obj$LSI_cluster = run_LSI(mtx, k = k)
}


if(cluster_method == 'scABC'){
  seurat.obj$scABC_cluster = run_scABC(mtx, k = k)
}

if(cluster_method == 'chromVar'){
  obj = run_chromVAR(mtx, genomeName)
  sele.cisTopic <- selectModel(obj, 
                               keepBinaryMatrix = F, keepModels = F)
  cell_topic <- t(modelMatSelection(sele.cisTopic, 'cell', 'Probability'))
  seurat.obj$chromVar_cluster = generalCluster(cell_topic, method = 'hclust', 
                                                             k = k)
}

saveRDS(seurat.obj, file = paste0(output_dir, '/seurat_obj_withCluster.rds'))

