library(data.table)
library(magrittr)
library(Seurat)
library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(BiocParallel)
#library(JASPAR2016)
#library(cisTopic)
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


filterMat <- function(atac.mtx, minFrac_in_cell = 0.01, min_depth = 200, max_depth = 100000){
  depth.cell = Matrix::colSums(atac.mtx)
  atac.mtx = atac.mtx[, depth.cell > min_depth & depth.cell < max_depth]
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
doBasicSeurat <- function(mtx, npc = 100, top.variable = 0.2, norm_by = 'log', 
                          doScale = T, doCenter = T, assay = 'ATAC', 
                          reg.var = 'nCount_ATAC'){
  
  # top.variabl -- use top most variable features
  if(doLog) mtx = round(log1p(mtx)/log(2))
  seurat.obj = CreateSeuratObject(mtx, project = 'scATAC', assay = assay,
                                  names.delim = '-')
  
  if(norm_by == 'log') seurat.obj@assays$ATAC@data <- log1p(seurat.obj@assays$ATAC@data)/log2
  if(norm_by == 'tf-idf') seurat.obj@assays$ATAC@data <- TF.IDF(seurat.obj@assays$ATAC@data)
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

regress_on_pca <- function(seurat.obj, reg.var = 'nCount_ATAC'){

  pcs = seurat.obj@reductions$pca@cell.embeddings
  pcs.reg = pcs
  for(i in 1:length(reg.var)){

    reg.var0 = seurat.obj[[reg.var[i]]][[1]]
    pcs.reg = apply(pcs.reg, 2, function(x) lm(x ~ reg.var0)$residual )

  }
   colnames(pcs.reg) = colnames(pcs)
  seurat.obj@reductions$pca@cell.embeddings = pcs.reg
  return(seurat.obj)
}


#could be normalized by log, tf-idf or none
doBasicSeurat_new <- function(mtx, npc = 50, top.variable = 0.2, 
                          doScale = T, doCenter = T, assay = 'ATAC',
                          reg.var = NULL, norm_by = 'log', project = 'scATAC'){

  # top.variabl -- use top most variable features
  seurat.obj = CreateSeuratObject(mtx, project = project, assay = assay,
                                  names.delim = '-')
 
  if(norm_by == 'log') seurat.obj@assays$ATAC@data <- log1p(seurat.obj@assays$ATAC@data)/log(2)
  if(norm_by == 'tf-idf') seurat.obj@assays$ATAC@data <- TF.IDF(seurat.obj@assays$ATAC@data)

  seurat.obj <- FindVariableFeatures(object = seurat.obj,
                                     selection.method = 'vst',
                                     nfeatures = floor(nrow(mtx) * top.variable))
  seurat.obj <- ScaleData(object = seurat.obj,
                          features = VariableFeatures(seurat.obj),
                          vars.to.regress = NULL, do.scale = doScale,
                          do.center = doCenter)


  seurat.obj <- RunPCA(object = seurat.obj,
                       features = VariableFeatures(object = seurat.obj),
                       verbose = FALSE, seed.use = 10, npc = npc)
  if(length(reg.var) > 0 ) seurat.obj = regress_on_pca(seurat.obj, reg.var)

 # seurat.obj <- RunLSI(seurat.obj, n = npc,
 #                      features = VariableFeatures(object = seurat.obj))

  return(seurat.obj)
}
doBasicSeurat_new = cmpfun(doBasicSeurat_new)

# assign gene to nearest peak and mark a gene if its tss within the peak
assignGene2Peak <- function(mtx, tss_ann){
  #tss_ann[, 'tss' := ifelse(strand == '+', start, end)]
  
  tss_ann[, 'tss' := start]
  peaks = tidyr::separate(data.table(x=rownames(mtx)),
                          col = x,
                          into = c('chr', 'start', 'end'))
  
  
  peaks$peak_name = rownames(mtx)
  class(peaks$start) = 'integer'
  class(peaks$end) = 'integer'
  
  
  chrs = unique(peaks$chr)
  peaks_ann = NULL
  for(chr0 in chrs){
    peaks0 = peaks[chr == chr0]
    genes0 = tss_ann[chr == chr0]
    peaks0[, 'id' := which.min(abs(genes0$tss - start/2 - end/2)), by = 'peak_name']
    peaks0[, 'gene_name' := genes0[id, gene_name]]
    peaks0$tss_name = ''
    for(i in 1:nrow(peaks0)){
      tss0 = genes0[tss <= peaks0$end[i] & tss >= peaks0$start[i]]
      if(nrow(tss0) > 0 ) {
        if(peaks0$gene_name[i] %in% tss0$gene_name) peaks0$gene_name[i] <- ''
        peaks0$tss_name[i] = paste(paste0(unique(tss0$gene_name), '-Tss'), 
                                                collapse = ',')
      }
    }
    
    peaks_ann = rbind(peaks_ann, peaks0)
  }
  peaks_ann[, 'id':= NULL]
  
  peaks_ann[, 'peak_new_name' := ifelse(!is.na(gene_name) & nchar(gene_name) > 1, 
                                        paste0(peak_name, ',', gene_name), peak_name)]
  
  
  peaks_ann[, 'peak_new_name' := ifelse(!is.na(tss_name) & nchar(tss_name) > 1, 
                                         paste0(peak_new_name, ',', tss_name), peak_new_name)]
  setkey(peaks_ann, peak_name)
  
  
  
  rownames(mtx) = peaks_ann[rownames(mtx)]$peak_new_name
  
  return(mtx)
  
  
}


# get the right resolution parater given number of expected cluster
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

## implement of generic clustering methods 
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

doDimReduction4mat <- function(mtx, max_pc = 20, doTSNE = F){
  ## DO SVD on tf_idf normalized matrix                             
  set.seed(1234)
  SVDtsne = irlba(mtx, max_pc, max_pc, maxit=1000)
  d_diagtsne = matrix(0, nrow=length(SVDtsne$d), ncol=length(SVDtsne$d))
  diag(d_diagtsne) = SVDtsne$d
  SVDtsne_vd = t(d_diagtsne %*% t(SVDtsne$v))
  rownames(SVDtsne_vd) = colnames(mtx)
  colnames(SVDtsne_vd) = paste0('pca_', 1:ncol(SVDtsne_vd))
  
  ## Run TSNE to 2 dimensions
  tsne_coords = NULL
  if(doTSNE){
    tsnetfidf = Rtsne(SVDtsne_vd, pca=F, perplexity=30, max_iter=5000)
    
    tsne_coords = as.data.frame(tsnetfidf$Y)
    colnames(tsne_coords) = c('tsne_1', 'tsne_2')
    rownames(tsne_coords) = colnames(mtx)
  }
  
  pca_coords = SVDtsne_vd
  return(list('pca_coords' = pca_coords, 'tsne_coords' = tsne_coords))
}


# the imput mtx is already filterd
run_scABC <- function(mtx, k = 5){
  weights = apply(mtx, 2, mean)
  landmarks = computeLandmarks(mtx, weights = weights, nCluster = k)
  labels = assign2landmarks(mtx, landmarks)
  return(labels)
}

run_chromVAR <- function(mtx, genomeName = 'BSgenome.Hsapiens.UCSC.hg38',
                         ncore = 3){
  
  register(MulticoreParam(3))
  if(!require(genomeName, character.only = T)) BiocManager::install(genomeName)  
  peaks = data.table('x' = rownames(mtx))
  peaks = tidyr::separate(peaks, col = 'x', into = c('chr', 'start', 'end'), sep = '-')
  peaks = GenomicRanges::makeGRangesFromDataFrame(peaks)
  
  frag.counts = SummarizedExperiment(assay = list(counts = mtx),
                                     rowRanges = peaks)
  frag.counts <- addGCBias(frag.counts, genome = genomeName)
  #motifs <- getJasparMotifs()
  library(chromVARmotifs) ## cisbp motif
  #motifs = ifelse(grepl(genomeName, pattern = 'hg'), human_pwms_v2, mouse_pwms_v2) 
  if(grepl(genomeName, pattern = 'hg')){
    motifs = human_pwms_v2
  }else{
    motifs = mouse_pwms_v2
  }
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
  if(!require(cisTopic)){
    if(!require(RcisTarge)) devtools::install_github("aertslab/RcisTarget")
    if(!require(AUCell)) devtools::install_github("aertslab/AUCell") 
    devtools::install_github("aertslab/cisTopic")
  }
  library(cisTopic)
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


# Do DA/DE with one cluster vs the rest clusters
# clusters are the data frame with <barcode> <cluster>
do_DA <- function(mtx_score, clusters, test = 'wilcox', 
                  only.pos = T, fdr = 0.05, topn = 10){
  clusters$cluster = as.character(clusters$cluster)
  cls = unique(clusters$cluster)
  res = NULL
  features = rownames(mtx_score)
  for(cluster0 in cls){
    bc0 = clusters[cluster == cluster0]$barcode
    mtx1 = mtx_score[, colnames(mtx_score) %in% bc0]
    mtx2 = mtx_score[, !colnames(mtx_score) %in% bc0]
    mu1 = sapply(1:length(features), function(x) mean(mtx1[x, ]))
    mu2 = sapply(1:length(features), function(x) mean(mtx2[x, ]))
    
   
    pvs = sapply(1:length(features), function(x) wilcox.test(mtx1[x, ], mtx2[x, ])$p.value )
    pvs.adj = p.adjust(pvs, method = 'fdr')
    res0 = data.table('feature' = features, 'cluster' = cluster0,
                      'mean1' = mu1, 'mean2' = mu2,
                       'pv' = pvs, 'pv_adjust' = pvs.adj)
    
    if(only.pos) res0 = res0[mean1 > 0]
    
    res0 = res0[order(pv_adjust), ]
    res0 = res0[pv_adjust <= fdr]
    
    if(nrow(res0) > topn) res0 = res0[1:topn, ]
    res = rbind(res, res0)
  }
  return(res)
}


#fg_genes: vector of forground genes
#bg_genes: background genes
#type: BP, CC, kegg
do_GO <- function(fg_genes, bg_genes, type = "BP", qCutoff = 0.05,
                 organism = c("mmu",  "hsa")) {
  if(organism =="mmu") {
    orgdb <- "org.Mm.eg.db"
    fromType = "SYMBOL"
    if(!require("org.Mm.eg.db")) BiocManager::install("org.Mm.eg.db")
    
  } else if(organism == "hsa") {
    orgdb <- "org.Hs.eg.db"
    if(!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
    fromType = "SYMBOL"
    
  }
  
  bg.df <- bitr(bg_genes, fromType = fromType,
                toType = c("SYMBOL", "ENTREZID"),
                OrgDb = orgdb)
  gene.df <- bitr(fg_genes, fromType = fromType,
                  toType = c("SYMBOL", "ENTREZID"),
                  OrgDb = orgdb)
  
  if(type == "kegg") {
    kegg_gene <- gene.df$ENTREZID
    kegg_bg <- bg.df$ENTREZID
    
    enrich_list <- enrichKEGG(
      gene          = kegg_gene,
      universe      = kegg_bg,
      organism      = organism,
      pAdjustMethod = "BH",
      qvalueCutoff  = qCutoff)
    
    
  }else {
    enrich_list <- enrichGO(gene        = gene.df$ENTREZID,
                            universe      = bg.df$ENTREZID,
                            OrgDb         = orgdb,
                            ont           = type,
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = qCutoff,
                            readable      = TRUE)
  }
  
  return(enrich_list@result[enrich_list@result$qvalue <= qCutoff, ])
}



#integrate two seurat object (with the same rownames)
integrateTwoSeurat <- function(seurat1, seurat2, npc = 50, reg.var = NULL,
                               reduction = 'pca'){
  seurat2$obj = 'obj2'
  seurat1$obj = 'obj1'
  comb.list = list('Obj1' = seurat1, 'Obj2' = seurat2)
  comb.anchors <- FindIntegrationAnchors(comb.list, scale = F, dims = 1:npc,
                                         anchor.features = 
                                        union(VariableFeatures(seurat1), 
                                              VariableFeatures(seurat2) ))
  comb.integrated <- IntegrateData(anchorset = comb.anchors, dims = 1:npc)
  
  # run standard workflow for integrated object
  DefaultAssay(object = comb.integrated) <- "integrated"
  comb.integrated <- ScaleData(object = comb.integrated, verbose = FALSE,
                               vars.to.regress = reg.var)
  if(reduction == 'pca'){
    comb.integrated <- RunPCA(object = comb.integrated, npcs = npc, verbose = FALSE)
  }
  if(reduction == 'lsi'){
    comb.integrated <- RunLSI(object = comb.integrated, n = npc, verbose = FALSE)
  }
  
  #comb.integrated <- RunTSNE(object = comb.integrated, reduction = reduction, 
  #                           dims = 1:npc)
  comb.integrated <- RunUMAP(object = comb.integrated, reduction = reduction, 
                             dims = 1:npc, verbose = F)
  return(comb.integrated)
}

#integrate seurat objects (with the same rownames)
integrateSeurats_atac <- function(seurat_list, npc = 50, reg.var = NULL,
                               reduction = 'pca'){
  #seurat2$sample = seurat2@project.name
  #seurat1$sample = seurat1@project.name
 
  ufeatures = lapply(seurat_list, function(x) VariableFeatures(x))
  ufeatures = unique(do.call('c', ufeatures))
  comb.anchors <- FindIntegrationAnchors(seurat_list, scale = F, dims = 1:npc,
                                         anchor.features = 6000)
  
  comb.integrated <- IntegrateData(anchorset = comb.anchors, dims = 1:npc)
  
  # run standard workflow for integrated object
  DefaultAssay(object = comb.integrated) <- "integrated"
  comb.integrated <- ScaleData(object = comb.integrated, verbose = FALSE,
                               vars.to.regress = NULL)
  if(reduction == 'pca'){
    comb.integrated <- RunPCA(object = comb.integrated, npcs = npc, verbose = FALSE)
    if(!is.null(reg.var)) comb.integrated <- regress_on_pca(comb.integrated, reg.var = reg.var)
  }
  if(reduction == 'lsi'){
    comb.integrated <- RunLSI(object = comb.integrated, n = npc, verbose = FALSE)
  }
  
  #comb.integrated <- RunTSNE(object = comb.integrated, reduction = reduction, 
  #                           dims = 1:npc)
  comb.integrated <- RunUMAP(object = comb.integrated, reduction = reduction, 
                             dims = 1:npc, verbose = F)
  return(comb.integrated)
}

# mtx: gene by cell matrix, or peak/bin by cell matrix
# seurat.obj: optional, if not provided, will construct one (which will clustering on pca1:50),
# if seurat.obj is provided, assume it did pca
# norm_mtx: normalized matrix, equals to log(mtx +1) if not provided
# extraDims: do extra dimension reduction on 10, 30, 100  etc
# subSamples: provided a subsample version
# vFeatures: given variable features; if NULL using default variable features
prepInput4Cello <- function(mtx, seurat.obj, norm_mtx = NULL,
                            cello.name = 'scATAC', assay = 'ATAC', 
                            resl4clust = 0.6, top.variable = 0.2, 
                            extraDims = c(10, 20, 30, 50, 80, 100),
                            subSamples = NULL, subCelloName = 'sub',
                            vars.to.regOnPca = NULL,
                            vFeatures = NULL){
  
  if(!is.null(vFeatures)) seurat.obj <- ScaleData(seurat.obj, features = vFeatures)
  if(is.null(vFeatures)) vFeatures = VariableFeatures(seurat.obj)  
  DefaultAssay(seurat.obj) = assay
  
  ndefault = ncol(seurat.obj@reductions$pca@cell.embeddings)
  
  seurat.obj <- RunPCA(seurat.obj, npcs = ndefault, verbose = F,
                       assay = assay, features = vFeatures)
  
  if(!is.null(vars.to.regOnPca)) seurat.obj = regress_on_pca(seurat.obj, vars.to.regOnPca)
  
  seurat.obj <- RunTSNE(seurat.obj, dims = 1:ndefault, check_duplicates = FALSE,
                        assay = assay)
  seurat.obj <- RunUMAP(seurat.obj, dims = 1:ndefault, verbose = F,
                        assay = assay)
  my_tsne_proj <- seurat.obj@reductions$tsne@cell.embeddings
  my_umap_proj <- seurat.obj@reductions$umap@cell.embeddings 
  
  
  
  if(ndefault < 100){
    
    seurat.obj <- RunPCA(seurat.obj, npcs = 100, verbose = F,
                         assay = assay, features = vFeatures)
    
    if(!is.null(vars.to.regOnPca)) seurat.obj = regress_on_pca(seurat.obj, vars.to.regOnPca)
    
    seurat.obj <- RunTSNE(seurat.obj, dims = 1:100, check_duplicates = FALSE,
                          assay = assay)
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:100, verbose = F,
                          assay = assay)
    my_tsne_proj100 <- seurat.obj@reductions$tsne@cell.embeddings
    my_umap_proj100 <- seurat.obj@reductions$umap@cell.embeddings 
    
    
    
  }
  
  
  
  #mtx = mtx[rownames(mtx) %in% rownames(seurat.obj), ]
  mtx = mtx[, colnames(mtx) %in% colnames(seurat.obj)]
  
  meta = seurat.obj@meta.data
  fmeta <- data.frame(symbol = rownames(mtx)) 
  rownames(fmeta) <- fmeta$symbol
  
  # preparing expression matrix
  if(is.null(norm_mtx)) norm_mtx = log2(mtx + 1)
  eset <- new("ExpressionSet",
              assayData = assayDataNew("environment", 
                                       exprs = mtx, 
                                       norm_exprs = norm_mtx),
              phenoData =  new("AnnotatedDataFrame", data = meta),
              featureData = new("AnnotatedDataFrame", data = fmeta))
  
  # Creating a cello 
  Cello <- setClass("Cello",
                    slots = c(
                      name = "character", # The name of cvis
                      idx = "numeric", # The index of global cds object
                      proj = "list", # The projections as a list of data frame
                      pmeta = "data.frame", # The local meta data
                      notes = "character" # Other information to display to the user
                    )
  )
  
  cello1 <- new("Cello", name = cello.name, idx = 1:ncol(mtx)) 
  
  my_pca_proj100 <- seurat.obj@reductions$pca@cell.embeddings
  my_tsne_proj100 <- seurat.obj@reductions$tsne@cell.embeddings
  my_umap_proj100 <- seurat.obj@reductions$umap@cell.embeddings 
  
  cello1@proj <- list("PCA" = my_pca_proj100) 
  if(is.null(extraDims)) extraDims = unique(c(ndefault, 100))
  
  for(dim0 in extraDims){
    if(dim0 != 100 & dim0 != ndefault){
      
      seurat.obj <- RunTSNE(seurat.obj, dims = 1:dim0, check_duplicates = FALSE,
                            assay = assay)
      seurat.obj <- RunUMAP(seurat.obj, dims = 1:dim0, verbose = F,
                            assay = assay)
      my_tsne_proj0 <- seurat.obj@reductions$tsne@cell.embeddings
      my_umap_proj0 <- seurat.obj@reductions$umap@cell.embeddings 
    }
    
    if(dim0 == ndefault){
      my_tsne_proj0 <- my_tsne_proj
      my_umap_proj0 <- my_umap_proj
    }
    
    if(dim0 == 100 & dim0 != ndefault){
      my_tsne_proj0 <- my_tsne_proj100
      my_umap_proj0 <- my_umap_proj100
    }
    
    cello1@proj[[paste0('t-SEN', dim0)]] <- my_tsne_proj0
    cello1@proj[[paste0('UMAP', dim0)]] <- my_umap_proj0
    
  }
  
  #cello1@proj = cello1@proj[order(names(cello1@proj))]
  
  clist <- list()
  clist[[cello.name]] <- cello1
  
  if(!is.null(subSamples)){
    ids = which(colnames(seurat.obj) %in% subSamples)
    cello2 <- new("Cello", name = cello.name, idx = ids) 
    
    
    cello2@proj <- list("PCA" = my_pca_proj[ids, ], 
                        't-SNE50' = my_tsne_proj50[ids, ],
                        "UMAP50" = my_umap_proj50[ids, ]) 
    
    cello.name.sub = paste0(cello.name, '_', subCelloName)
    clist[[cello.name.sub]] <- cello2
  }
  
  
  return(list('eset' = eset, 'clist' = clist))
}

