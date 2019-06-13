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
genome_name = args[5]

mtx = read_mtx_scATACpro(mtx_file)
mtx = filterMat(mtx)


if(file.exists(paste0(output_dir, '/seurat_obj_withCluster.rds'))){
  seurat.obj = readRDS(paste0(output_dir, '/seurat_obj_withCluster.rds'))
}else if(file.exists(paste0(output_dir, '/seurat_obj_withDR.rds'))){
  seurat.obj = readRDS(paste0(output_dir, '/seurat_obj_withDR.rds'))
}else{
  seurat.obj = doBasicSeurat(mtx, npc = 100, top.variable = 0.2)
  
  seurat.obj = RunTSNE(seurat.obj, dims = 1:50)
  seurat.obj = RunUMAP(seurat.obj, dims = 1:50, verbose = F)
  saveRDS(seurat.obj, file = paste0(output_dir, '/seurat_obj_withDR.rds'))
}




## clustering
if(cluster_method == 'seurat'){
  ## seurat implemented louvain algorithm
  seurat.obj = FindNeighbors(seurat.obj, reduction = 'pca')
  resl = queryResolution4Seurat(seurat.obj, reduction = 'pca', npc = 50, k = k,
                                min_resl = 0.05)
  seurat.obj = FindClusters(seurat.obj, resolution = resl)
  seurat.obj$seurat_clusters = seurat.obj@active.ident
  seurat.obj$active_clusters = seurat.obj$seurat_clusters
}

if(cluster_method == 'cisTopic'){
  nc = detectCores()
  cis.obj = run_cisTopic(mtx, nCores = max(1, nc - 1))
  sele.cisTopic <- selectModel(cis.obj, 
                               keepBinaryMatrix = F, keepModels = F)
  cell_topic <- t(modelMatSelection(sele.cisTopic, 'cell', 'Probability'))
  seurat.obj$cisTopic_clusters = generalCluster(cell_topic, method = 'hclust', k = k)
  seurat.obj$active_clusters = seurat.obj$cisTopic_clusters
}

if(cluster_method == 'LSI'){
  seurat.obj$LSI_clusters = run_LSI(mtx, k = k)
  seurat.obj$active_clusters = seurat.obj$LSI_clusters
}


if(cluster_method == 'scABC'){
  seurat.obj$scABC_clusters = run_scABC(mtx, k = k)
  seurat.obj$active_clusters = seurat.obj$scABC_clusters
}

if(cluster_method == 'chromVar'){
  genomeName = 'BSgenome.Hsapiens.UCSC.hg38'
  if(grepl(genome_name, pattern = '38'))genomeName = 'BSgenome.Hsapiens.UCSC.hg38'
  if(grepl(genome_name, pattern = '19'))genomeName = 'BSgenome.Hsapiens.UCSC.hg19'
  if(grepl(genome_name, pattern = 'mm9'))genomeName = 'BSgenome.Mmusculus.UCSC.mm9'
  if(grepl(genome_name, pattern = 'mm10'))genomeName = 'BSgenome.Mmusculus.UCSC.mm10'
  nc = detectCores()
  obj = run_chromVAR(mtx, genomeName, max(1, nc - 1))
  saveRDS(obj, file = paste0(output_dir, '/chromVar_obj.rds'))
  pca_coords = doDimReduction4mat(obj@assays$data$z)[[1]]
  
  seurat.obj$chromVar_clusters = cutree(hclust(dist(pca_coords)), k = k)
  seurat.obj$active_clusters = seurat.obj$chromVar_clusters
  
  
}

saveRDS(seurat.obj, file = paste0(output_dir, '/seurat_obj_withCluster.rds'))



