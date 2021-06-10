source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local('dsAnalysis_utilities.R')

library(ggplot2)
library(parallel)
args = commandArgs(T)
mtx_file = args[1]
cluster_method = args[2]
k = args[3]
output_dir = args[4]
genome_name = args[5]
tss_path = args[6]
norm_by = args[7]
REDUCTION = args[8]
nREDUCTION = as.integer(args[9])
top_variable_features = as.numeric(args[10])
qc_stat_file = args[11]

if(grepl(mtx_file, pattern = '.rds', fix = T)) {
    mtx = readRDS(mtx_file)
    if(any(class(mtx) == 'Seurat')) mtx = mtx@assays$ATAC@counts
    rnames = rownames(mtx)
    new.rnames = sapply(rnames, function(x) gsub('_', '-', x))
    names(new.rnames) = NULL
    rownames(mtx) <- new.rnames
}else{
    mtx = read_mtx_scATACpro(mtx_file)
}

tss_ann <- fread(tss_path, header = F)
#names(tss_ann)[c(1:4,7)] <- c('chr', 'start', 'end', 'gene_name', 'gene_type')
#tss_ann <- tss_ann[gene_type %in% c('miRNA', 'lincRNA', 'protein_coding'), ]
names(tss_ann)[c(1:4)] <- c('chr', 'start', 'end', 'gene_name')
mtx = assignGene2Peak(mtx, tss_ann)


## remove peaks than are less freqent than 0.5% of cells
rfreqs = Matrix::rowMeans(mtx > 0)
mtx = mtx[rfreqs > 0.005, ]
cfreqs = Matrix::colMeans(mtx > 0)
mtx = mtx[, cfreqs > 0]

seurat.obj = runSeurat_Atac(mtx, npc = nREDUCTION, norm_by = norm_by, 
                               top_variable_features = top_variable_features, reg.var = 'nCount_ATAC')

## add qc stat to each cell
if(file.exists(qc_stat_file)){
    qc_singlecell = fread(qc_stat_file)
    qc_singlecell = qc_singlecell[bc %in% colnames(seurat.obj)]
    qc_singlecell = data.frame(qc_singlecell)
    rownames(qc_singlecell) = qc_singlecell$bc
    qc_singlecell$bc = NULL
    names(qc_singlecell) =  c("total.unique.frags", "frac.mito",  "frac.peak",
                             "frac.promoter", "frac.tss", "frac.enhancer", "tss_enrich_score")
    seurat.obj <- AddMetaData(seurat.obj, metadata = qc_singlecell)
}




if(REDUCTION != 'lda'){
    seurat.obj = RunTSNE(seurat.obj, dims = 1:nREDUCTION, reduction = 'pca', check_duplicates = FALSE)
    seurat.obj = RunUMAP(seurat.obj, dims = 1:nREDUCTION, reduction = 'pca', verbose = F)
}

saveRDS(seurat.obj, file = paste0(output_dir, '/seurat_obj.rds'))




## clustering
if(cluster_method == 'seurat'){
  ## seurat implemented louvain algorithm
  seurat.obj = FindNeighbors(seurat.obj, reduction = 'pca', dims = 1:30, k.param = 20)
  if (toupper(k) == 'NULL' || k == '0'){
    resl = 0.2
  }else{
    k = as.numeric(k)
    ki = as.integer(k)
    if(k >= 2 & k == ki) {
        resl = queryResolution4Seurat(seurat.obj, reduction = 'pca', npc = nREDUCTION, k = ki,
                            min_resl = 0.01)
    }else{
        resl = k
    }
  }
  seurat.obj = FindClusters(seurat.obj, resolution = resl)
  seurat.obj$active_clusters = seurat.obj$seurat_clusters
}

k = as.integer(k)

if(grepl(REDUCTION, pattern = 'lda', ignore.case = T)){
    
  cis.obj = run_cisTopic(mtx, nCores = 2, topic = nREDUCTION)
  sele.cisTopic <- cisTopic::selectModel(cis.obj, select = nREDUCTION, 
                               keepBinaryMatrix = F, keepModels = F)
  cell_topic <- t(modelMatSelection(sele.cisTopic, 'cell', 'Probability'))
  seurat.obj[['lda']] <- CreateDimReducObject(embeddings = cell_topic, 
                                              key = 'Topic', assay = DefaultAssay(seurat.obj))
    seurat.obj = RunTSNE(seurat.obj, dims = 1:ncol(cell_topic), reduction = 'lda', check_duplicates = FALSE)
    seurat.obj = RunUMAP(seurat.obj, dims = 1:ncol(cell_topic), reduction = 'lda', verbose = F)
}

if(cluster_method == 'cisTopic' ){
  nc = detectCores()
  topic = unique(c(10, 20, 30, 50, 80, 100, nREDUCTION))
  cis.obj = run_cisTopic(mtx, nCores = min(10, nc), topic = topic)
  sele.cisTopic <- cisTopic::selectModel(cis.obj,
                               keepBinaryMatrix = F, keepModels = F)
  cell_topic <- t(modelMatSelection(sele.cisTopic, 'cell', 'Probability'))
  seurat.obj$cisTopic_clusters = generalCluster(cell_topic, method = 'hclust', k = k)
  seurat.obj$active_clusters = seurat.obj$cisTopic_clusters
  #seurat.obj[['lda']] <- CreateDimReducObject(embeddings = cell_topic, 
  #                                            key = 'Topic', assay = DefaultAssay(seurat.obj))
  #  seurat.obj = RunTSNE(seurat.obj, dims = 1:ncol(cell_topic), reduction = 'lda', check_duplicates = FALSE)
  #  seurat.obj = RunUMAP(seurat.obj, dims = 1:ncol(cell_topic), reduction = 'lda', verbose = F)
  #saveRDS(cell_topic, file = paste0(output_dir, '/cell_topic_obj.rds'))
}

if(cluster_method == 'LSI'){
  seurat.obj$LSA_clusters = run_LSI(mtx, k = k)
  seurat.obj$active_clusters = seurat.obj$LSI_clusters
}

if(cluster_method == 'SCRAT'){
  seurat.obj$SCRAT_clusters = run_scrat(mtx, k = k)
  seurat.obj$active_clusters = seurat.obj$SCRAT_clusters
}

if(cluster_method == 'kmeans'){
  seurat.obj$kmeans_clusters = generalCluster(seurat.obj@reductions$pca@cell.embeddings, k = k)
  seurat.obj$active_clusters = seurat.obj$kmeans_clusters
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
  seurat.obj$active_clusters = as.factor(seurat.obj$chromVar_clusters)
  
}

saveRDS(seurat.obj, file = paste0(output_dir, '/seurat_obj.rds'))

#output bacrcode cluster information
bc_cls = data.table('Barcode' = rownames(seurat.obj@meta.data), 'Cluster' = seurat.obj@meta.data$active_clusters)
setkey(bc_cls, Cluster)

write.table(bc_cls, file = paste0(output_dir, '/cell_cluster_table.tsv'), sep = '\t',
            quote = F, row.names = F)

cg <- DimPlot(seurat.obj, reduction = 'umap', group.by = 'active_clusters', label = T) + theme(legend.text = element_text(size = 17))
    
   pfname = paste0(output_dir, '/umap_clusters.eps')
   ggsave(cg, file = pfname, device = 'eps', width = 6, height = 6)

