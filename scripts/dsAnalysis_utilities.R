library(data.table)
library(magrittr)
library(Seurat)
library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(BiocParallel)
library(JASPAR2016)
library(cisTopic)
library(compiler)
library(readr)

## do reverse complemente of a DNA sequence
rev.comp <- function(x, rev=TRUE){
  x<-toupper(x)
  y<-rep("N",nchar(x))
  xx<-unlist(strsplit(x,NULL))
  for (bbb in 1:nchar(x))
  {
    if(xx[bbb]=="A") y[bbb]<-"T"    
    if(xx[bbb]=="C") y[bbb]<-"G"    
    if(xx[bbb]=="G") y[bbb]<-"C"    
    if(xx[bbb]=="T") y[bbb]<-"A"
  }
  if(rev==FALSE) 
  {
    for(ccc in (1:nchar(x)))
    {
      if(ccc==1) yy<-y[ccc] else yy<-paste(yy,y[ccc],sep="")
    }
  }
  if(rev==T)
  {
    zz<-rep(NA,nchar(x))
    for(ccc in (1:nchar(x)))
    {
      zz[ccc]<-y[nchar(x)+1-ccc]
      if(ccc==1) yy<-zz[ccc] else yy<-paste(yy,zz[ccc],sep="")
    }
  }
  return(yy)  
}


read_mtx_scATACpro <- function(mtx_path){
  #mtx_path <- paste0(dirt, "matrix.mtx")
  mtx.dir = dirname(mtx_path)
  feature_path <- paste0(mtx.dir, "/features.txt")
  barcode_path <- paste0(mtx.dir, "/barcodes.txt")
  
  
  features <- fread(feature_path, header = F)
  barcodes <- fread(barcode_path, header = F)
  
  mtx <-  Matrix::readMM(mtx_path) %>%
    magrittr::set_rownames(features$V1)%>%
    magrittr::set_colnames(barcodes$V1) 
  
  return(mtx)
}


filterMat <- function(atac.mtx, minFrac_in_cell = 0.01, min_depth = 1000){
  depth.cell = Matrix::colSums(atac.mtx)
  atac.mtx = atac.mtx[, depth.cell > min_depth]
  frac.in.cell = Matrix::rowSums(atac.mtx > 0)
  atac.mtx = atac.mtx[frac.in.cell > minFrac_in_cell, ]
  return(atac.mtx)
}

## tf-idf normalization
atac_tfidf = function(atac_matrix, site_frequency_threshold=0.03) {
  num_cells_ncounted = Matrix::rowSums(atac_matrix)
  threshold = ncol(atac_matrix) * site_frequency_threshold
  
  ncounts = atac_matrix[num_cells_ncounted >= threshold,]
  
  ## Normalize the data with TF-IDF
  nfreqs = t(t(ncounts) / Matrix::colSums(ncounts))
  tf_idf_counts = nfreqs * log(1 + ncol(ncounts) / Matrix::rowSums(ncounts))
  
  return(list('tfidf_matrix' = tf_idf_counts, 'filtered_matrix' = ncounts))
}


# do normalization, pca using Seurat
doBasicSeurat <- function(mtx, npc = 100, top.variable = 0.2, doLog = T, 
                          doScale = T, doCenter = T, assay = 'ATAC', 
                          reg.var = 'nCount_ATAC'){
  
  # top.variabl -- use top most variable features
  if(doLog) mtx = log2(1 + mtx)
  seurat.obj = CreateSeuratObject(mtx, project = 'scATAC', assay = assay,
                                  names.delim = '-')
  cell.names = colnames(mtx)
  
  #seurat.obj <- NormalizeData(seurat.obj, normalization.method = 'LogNormalize',
  #                            scale.factor = 1e4)
  
  seurat.obj <- FindVariableFeatures(object = seurat.obj, 
                                     selection.method = 'vst', 
                                     nfeatures = floor(nrow(mtx) * top.variable))
  seurat.obj <- ScaleData(object = seurat.obj, 
                          features = VariableFeatures(seurat.obj), 
                          vars.to.regress = reg.var, do.scale = doScale, do.center = doCenter)
  
  
  seurat.obj <- RunPCA(object = seurat.obj, 
                       features = VariableFeatures(object = seurat.obj),
                       verbose = FALSE, seed.use = 10, npc = npc)
  
  return(seurat.obj)
}
doBasicSeurat = cmpfun(doBasicSeurat)


queryResolution4Seurat <- function(seurat.obj, k = 10, reduction = 'umap', npc = 20, 
                                   min_resl = 0.1, max_resl = 1, max_iter = 15, doPCA = F){
  max.dim = ifelse(reduction == 'pca', npc, 2)
  if(doPCA) {
    seurat.obj <- RunPCA(seurat.obj, npcs = npc, verbose = F)
    seurat.obj <- RunTSNE(seurat.obj, dims = 1:npc, verbose = F)
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:npc, verbose = F)
  }
  
  
  seurat.obj <- FindNeighbors(seurat.obj, reduction = reduction, verbose = F, dims = 1:max.dim)
  tmp.cluster1 <- FindClusters(seurat.obj, resolution = min_resl)@active.ident
  tmp.cluster2 <- FindClusters(seurat.obj, resolution = max_resl)@active.ident
  
  
  
  
  len1 = length(levels(tmp.cluster1))
  len2 = length(levels(tmp.cluster2))
  
  k1 = k2 = 0
  while(len1 > k ){
    
    k1 = k1 + 1
    message('min_resl too large, trying to divided it by  2')
    min_resl = min_resl/2
    tmp.cluster1 <- FindClusters(seurat.obj, resolution = min_resl)@active.ident
    len1 = length(levels(tmp.cluster1))
    if(k1 == 10) stop('Please specify a much smaller min_res')
  }
  
  while(len2 < k){
    k2 = k2 + 1
    message('max_resl too small, trying to multiply it by 2')
    max_resl = max_resl * 2
    tmp.cluster2 <- FindClusters(seurat.obj, resolution = max_resl)@active.ident
    len2 = length(levels(tmp.cluster2))
    if(k2 == 10) stop('Please specify a much bigger max_res')
  }
  if(len1 == k) {
    return(min_resl)
  }
  
  if(len2 == k) {
    return(max_resl)
  }
  
  # repeat in other case
  
  i = 0
  repeat{
    i = i + 1
    resl0 = min_resl/2 + max_resl/2
    
    tmp.cluster <- FindClusters(seurat.obj, resolution = resl0)@active.ident
    
    len = length(levels(tmp.cluster)) 
    if(len == k){
      return(resl0)
    }
    if(len < k){
      min_resl = resl0
      len1 = len
    }
    if(len > k){
      max_resl = resl0
      len2 = len
    }
    if(i == max_iter) break
  }
  return(resl0)
}
queryResolution4Seurat = cmpfun(queryResolution4Seurat)

## implement of generic clustering methods ####
generalCluster <- function(reduced.mtx, method = 'hclust', k = 5){
  if(method == 'kmeans'){
    res = kmeans(reduced.mtx, centers = k)
    cl.label = res$cluster
  }
  
  if(method == 'hclust'){
    if(is.null(k)) stop('Need specify k: the number of cluster')
    d <- dist(reduced.mtx, method = "euclidean") # distance matrix
    fit <- hclust(d, method = "ward.D")
    cl.label<- cutree(fit, k = k) # cut tree into k clusters
  }
  
  if(method == 'mclust'){
    
    fit <- Mclust(reduced.mtx, G = k)
    cl.label = fit$classification
  }
  return(cl.label)
}

# the imput mtx is already filterd
run_scABC <- function(mtx, k = 5){
  weights = apply(mtx, 2, mean)
  landmarks = computeLandmarks(mtx, weights = weights, nCluster = k)
  labels = assign2landmarks(mtx, landmarks)
  return(labels)
}

run_chromVAR <- function(mtx, genomeName = 'BSgenome.Hsapiens.UCSC.hg19'){
  
  register(MulticoreParam(4))
  
  peaks = data.table('x' = rownames(mtx))
  peaks = tidyr::separate(peaks, col = 'x', into = c('chr', 'start', 'end'))
  peaks = GenomicRanges::makeGRangesFromDataFrame(peaks)
  
  frag.counts = SummarizedExperiment(assay = list(counts = mtx),
                                     rowRanges = peaks)
  frag.counts <- addGCBias(frag.counts, genome = genomeName)
  motifs <- getJasparMotifs()
  motif_ix <- matchMotifs(motifs, frag.counts,
                          genome = genomeName)
  dev <- computeDeviations(object = frag.counts, 
                           annotations = motif_ix)
  bg <- getBackgroundPeaks(object = frag.counts)
  
  dev <- computeDeviations(object = frag.counts, annotations = motif_ix,
                           background_peaks = bg)
  
  #motif.zscore = dev@assays$data$z
  return(dev)
}

run_LSI <- function(mtx, ncell.peak = 150,  max_pc = 10, k = 5){
  mtx = (mtx > 0)
  mtx = 1 * mtx
  
  num_cells.peak = Matrix::rowSums(mtx)
  ncounts = mtx[num_cells.peak >= ncell.peak,]
  
  ## Normalize the data with TF-IDF
  nfreqs = t(t(ncounts) / Matrix::colSums(ncounts))
  tf_idf_counts = nfreqs * log(1 + ncol(ncounts) / Matrix::rowSums(ncounts))
  
  ## DO SVD on tf_idf normalized matrix                             
  set.seed(0)
  SVDtsne = irlba(tf_idf_counts, max_pc, max_pc, maxit=1000)
  d_diagtsne = matrix(0, nrow=length(SVDtsne$d), ncol=length(SVDtsne$d))
  diag(d_diagtsne) = SVDtsne$d
  SVDtsne_vd = t(d_diagtsne %*% t(SVDtsne$v))
  rownames(SVDtsne_vd) = colnames(mtx)
  colnames(SVDtsne_vd) = paste0('pca_', 1:ncol(SVDtsne_vd))
  
  ## Run TSNE to 2 dimensions
  if(F){
    tsnetfidf = Rtsne(SVDtsne_vd, pca=F, perplexity=30, max_iter=5000)
    
    tsne_coords = as.data.frame(tsnetfidf$Y)
    colnames(tsne_coords) = c('tsne_1', 'tsne_2')
    rownames(tsne_coords) = colnames(ncounts) 
  }
  
  pca_coords = SVDtsne_vd[, 2:max_pc]
  cl.labels = generalCluster(pca_coords, k = k, method = 'hclust')
  return(cl.labels)
}

run_cisTopic <- function(mtx, nCores = 4){
  # prepare the right format of rownames
  rnames = data.table('region' = rownames(mtx))
  tmp = tidyr::separate(rnames, col = 'region', into = c('chr', 'start', 'end'))
  rnames = paste0(tmp$chr, ':', tmp$start, '-', tmp$end)
  rownames(mtx) = rnames
  
  cisTopicObject <- createcisTopicObject(mtx, project.name='scATAC')
  cisTopicObject <- runModels(cisTopicObject, topic = c(5, 10, 15, 20, 30, 40,
                                                        50, 60, 70, 80, 90, 100), seed = 987, nCores = nCores, 
                              burnin = 120, iterations = 150, addModels = T)
  #cisTopicObject <- selectModel(cisTopicObject, keepBinarymatrix = F, keepModels = F)
  #cellassign <- t(modelMatSelection(cisTopicObject, 'cell', 'Probability'))
  return(cisTopicObject)
}


