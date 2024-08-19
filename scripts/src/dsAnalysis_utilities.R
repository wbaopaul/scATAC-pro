library(data.table)
library(magrittr)
library(Seurat)
#library(chromVAR)
#library(motifmatchr)
library(SummarizedExperiment)
library(BiocParallel)
library(compiler)
library(readr)
library(matrixStats)
library(GenomicRanges)
library(edgeR)
library(mclust)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(Matrix)

TF_IDF <- function (data, verbose = T) 
{
    if (class(x = data) == "data.frame") {
        data <- as.matrix(x = data)
    }
    if (class(x = data) != "dgCMatrix") {
        data <- as(object = data, Class = "dgCMatrix")
    }
    if (verbose) {
        message("Performing TF_IDF normalization")
    }
    npeaks <- colSums(x = data) + 1
    tf <- t(x = t(x = data)/npeaks)
    idf <- ncol(x = data)/(rowSums(x = data) + 1)
    idf <- log(1 + idf)
    norm.data <- Diagonal(n = length(x = idf), x = idf) %*% tf
    #norm.data[which(x = is.na(x = norm.data))] <- 0
    rownames(norm.data) <- rownames(data)
    return(norm.data)
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

## cbind saprse mtx with the union of the rownames
cBind_union_features <- function(mat_list){
    ff0 = ff = rownames(mat_list[[1]])
    for(i in 2:length(mat_list)){
      ff = unique(union(ff, rownames(mat_list[[i]])))
    }
    if(length(ff) == length(ff0)) return(do.call('cbind', mat_list))
   ## make a mtx with full features
    mat_union = list()
    for(i in 1:length(mat_list)){
      mtx0 = mat_list[[i]]
      ff0 = setdiff(ff, rownames(mtx0))
      if(length(ff0) > 0 ) {
        tmp = as(matrix(0, length(ff0), ncol(mtx0)), "sparseMatrix")
        rownames(tmp) = ff0
        mtx0 = rbind(mtx0, tmp)
      }
      mat_union[[i]] = mtx0[order(rownames(mtx0)), ]
    }
    return(do.call('cbind', mat_union))
}

filterMat <- function(atac.mtx, minFrac_in_cell = 0.01, min_depth = 200, 
                      max_depth = 100000){
  depth.cell = Matrix::colSums(atac.mtx)
  atac.mtx = atac.mtx[, depth.cell > min_depth & depth.cell < max_depth]
  frac.in.cell = Matrix::rowMeans(atac.mtx > 0)
  atac.mtx = atac.mtx[frac.in.cell > minFrac_in_cell, ]
  depth.cell = Matrix::colSums(atac.mtx)
  atac.mtx = atac.mtx[, depth.cell > 0] ## remove cells without any read overlap filtered peaks
  return(atac.mtx)
}

## tf-idf normalization
atac_tfidf = function(atac_matrix, site_frequency_threshold=0.03) {
  num_cells_ncounted = Matrix::rowSums(atac_matrix)
  threshold = ncol(atac_matrix) * site_frequency_threshold
  
  ncounts = atac_matrix[num_cells_ncounted >= threshold,]
  
  ## Normalize the data with TF_IDF
  nfreqs = t(t(ncounts) / Matrix::colSums(ncounts))
  tf_idf_counts = nfreqs * log(1 + ncol(ncounts) / Matrix::rowSums(ncounts))
  
  return(list('tfidf_matrix' = tf_idf_counts, 'filtered_matrix' = ncounts))
}



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
runSeurat_Atac <- function(mtx, npc = 50, top_variable_features = 0.2, 
                          doScale = T, doCenter = T, assay = 'ATAC',
                          reg.var = NULL, norm_by = 'log', project = 'scATAC'){

  # top.variabl -- use top most variable features
  seurat.obj = CreateSeuratObject(mtx, project = project, assay = assay,
                                  names.delim = '-', min.cells = 1,
                                  min.features = 1)
  
  if(class(seurat.obj[[assay]]) == 'Assay5'){
    if(norm_by == 'log') seurat.obj[[assay]]$data <- log1p(seurat.obj[[assay]]$counts)/log(2)
    if(norm_by == 'tf-idf') seurat.obj[[assay]]$data <- TF_IDF(seurat.obj[[assay]]$counts)
  }else{
    if(norm_by == 'log') seurat.obj[[assay]]@data <- log1p(seurat.obj[[assay]]@counts)/log(2)
    if(norm_by == 'tf-idf') seurat.obj[[assay]]@data <- TF_IDF(seurat.obj[[assay]]@counts)
  }
  
  nvap = ifelse(top_variable_features > 1, top_variable_features, 
                floor(nrow(mtx) * top_variable_features))

  seurat.obj <- FindVariableFeatures(object = seurat.obj,
                                     selection.method = 'vst',
                                     nfeatures = nvap)
  
  ## remove variable features only accessible in less than 1% of cells
 
  if(class(seurat.obj[[assay]]) == 'Assay5'){
    mtx = seurat.obj[[assay]]$counts
  }else{
    mtx = seurat.obj[[assay]]@counts
  }
  rs = Matrix::rowMeans(mtx > 0)
  rare.features = names(which(rs < 0.01))
  vaps = VariableFeatures(seurat.obj)
  vaps = setdiff(vaps, rare.features)
  niter = 0
  while(length(vaps) < 500 & nvap > 500){
    niter = niter + 1
    nvap = nvap + 2000
    seurat.obj <- FindVariableFeatures(object = seurat.obj,
                                       selection.method = 'vst',
                                       nfeatures = min(nvap, nrow(seurat.obj)))
    vaps = VariableFeatures(seurat.obj)
    vaps = setdiff(vaps, rare.features)
    if(niter >= 5) break 
  }

  if(length(vaps) < 1000) {
    message('Most variable features were too rare and were filtered, \n
            I am gonna select from them top 1000 less rare!')
    vaps = VariableFeatures(seurat.obj)
    vaps = names(sort(rs[vaps], decreasing = T))[1:1000]
  }
  VariableFeatures(seurat.obj) <- vaps
  
  ## redo normalization using vap if norm by tf-idf
  if(norm_by == 'tf-idf'){
    mtx.norm = TF_IDF(mtx[vaps, ])
    tmp <- mtx[setdiff(rownames(mtx), vaps), ]
    data0 <- rbind(mtx.norm, tmp)
    if(class(seurat.obj[[assay]]) == 'Assay5'){
      seurat.obj[[assay]]$data = data0[rownames(mtx), ]
    }else{
      seurat.obj[[assay]]@data = data0[rownames(mtx), ]
    }
    
    
    rm(data0, tmp, mtx.norm)
  }

  seurat.obj <- ScaleData(object = seurat.obj,
                          features = VariableFeatures(seurat.obj),
                          vars.to.regress = NULL, do.scale = doScale,
                          do.center = doCenter)


  seurat.obj <- RunPCA(object = seurat.obj,
                       features = VariableFeatures(object = seurat.obj),
                       verbose = FALSE, seed.use = 10, npc = npc)
  if(length(reg.var) > 0 ) seurat.obj = regress_on_pca(seurat.obj, reg.var)

  if(norm_by == 'tf-idf'){
  ## redo normalization using all features
     mtx.norm = TF_IDF(mtx)
     if(class(seurat.obj[[assay]]) == 'Assay5'){
       seurat.obj[[assay]]$data = mtx.norm
     }else{
       seurat.obj[[assay]]@data = mtx.norm
     }
     
  }

  return(seurat.obj)
}
runSeurat_Atac = cmpfun(runSeurat_Atac)

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
    peaks0[, 'id' := which.min(abs(genes0$tss - start/2 - end/2)), by = peak_name]
    peaks0[, 'gene_name' := genes0[id, gene_name]]
    peaks0[, 'dist0' := min(abs(genes0$tss - start/2 -end/2)), by = peak_name]
    peaks0[, 'gene_name' := ifelse(dist0 > 100000, '', gene_name)]
    peaks0$tss_name = ''
    for(i in 1:nrow(peaks0)){
      tss0 = genes0[tss <= (peaks0$end[i] + 1000) & tss >= (peaks0$start[i] - 1000)]
      if(nrow(tss0) > 0 ) {
        if(peaks0$gene_name[i] %in% tss0$gene_name) peaks0$gene_name[i] <- ''
        peaks0$tss_name[i] = paste(paste0(unique(tss0$gene_name), '-Tss'), 
                                                collapse = ',')
      }
    }
    
    peaks_ann = rbind(peaks_ann, peaks0)
  }
  peaks_ann[, 'id':= NULL]
  peaks_ann[, 'dist0':= NULL]
  
  peaks_ann[, 'peak_new_name' := ifelse(!is.na(gene_name) & nchar(gene_name) > 1, 
                                        paste0(peak_name, ',', gene_name), peak_name)]
  
  
  peaks_ann[, 'peak_new_name' := ifelse(!is.na(tss_name) & nchar(tss_name) > 1, 
                                         paste0(peak_new_name, ',', tss_name), peak_new_name)]
  setkey(peaks_ann, peak_name)
  
  
  
  rownames(mtx) = peaks_ann[rownames(mtx)]$peak_new_name
  
  return(mtx)
  
  
}



# assign gene to nearest peak and mark a gene if its tss within the peak
# input peak_coords with chr-start-end, format
assignGene2Peak_coords <- function(peak_coords, tss_ann){
  tss_ann[, 'tss' := start]
  peaks = tidyr::separate(data.table(x = peak_coords),
                          col = x,
                          into = c('chr', 'start', 'end'))
  
  
  peaks$peak_name = peak_coords
  class(peaks$start) = 'integer'
  class(peaks$end) = 'integer'
  
  geneTssInPeak <- function(tss_ids, genes0){
    if(!is.na(tss))
      rr = genes0[tss <= end & tss >= start]$gene_name
    return(paste(rr, collapse = ','))
  }
  
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
  
  return(peaks_ann[peak_coords, ]$peak_new_name)
  
}


# get the right resolution parater given number of expected cluster
queryResolution4Seurat <- function(seurat.obj, k = 10, reduction = 'umap', npc = 20, 
                                   min_resl = 0.1, max_resl = 1, max_iter = 15, doPCA = F){
  max.dim = ifelse(reduction == 'pca', npc, 2)
  if(doPCA) {
    seurat.obj <- RunPCA(seurat.obj, npcs = npc, verbose = F)
    seurat.obj <- RunTSNE(seurat.obj, dims = 1:npc, verbose = F, check_duplicates = F)
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:npc, verbose = F)
  }
  
  
  seurat.obj <- FindNeighbors(seurat.obj, reduction = reduction, verbose = F, dims = 1:max.dim)
  tmp.cluster1 <- FindClusters(seurat.obj, resolution = min_resl, verbose = F)@active.ident
  tmp.cluster2 <- FindClusters(seurat.obj, resolution = max_resl, verbose = F)@active.ident
  
  
  
  
  len1 = length(levels(tmp.cluster1))
  len2 = length(levels(tmp.cluster2))
  
  k1 = k2 = 0
  while(len1 > k ){
    
    k1 = k1 + 1
    message('min_resl too large, trying to divided it by  2')
    min_resl = min_resl/2
    tmp.cluster1 <- FindClusters(seurat.obj, resolution = min_resl, verbose=F)@active.ident
    len1 = length(levels(tmp.cluster1))
    if(k1 == 10) stop('Please specify a much smaller min_res')
  }
  
  while(len2 < k){
    k2 = k2 + 1
    message('max_resl too small, trying to multiply it by 2')
    max_resl = max_resl * 2
    tmp.cluster2 <- FindClusters(seurat.obj, resolution = max_resl, verbose = F)@active.ident
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
    
    tmp.cluster <- FindClusters(seurat.obj, resolution = resl0, verbose = F)@active.ident
    
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
  library(chromVAR)
  library(motifmatchr)
  rs = Matrix::rowSums(mtx)
  mtx = mtx[rs > 0, ]  
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
  
  ## Normalize the data with TF_IDF
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

run_scrat <- function(mtx, reduction = 'pca', max_pc = 20, method = 'mclust', k = 10){
  # stadardized features per cell
  mtx = scale(mtx, center = T, scale = T)
  if(reduction == 'pca'){
    reduced.mtx = doDimReduction4mat(mtx, max_pc = max_pc)[[1]]
    reduced.mtx = reduced.mtx[, -1]
  }
  if(reduction == 'tsne'){
    reduced.mtx = doDimReduction4mat(mtx, max_pc = max_pc, doTSNE = T)[[2]]
  }
  
  cl.label = generalCluster(reduced.mtx, k = k, method = method)
  return(cl.label)
}


run_cisTopic <- function(mtx, nCores = 4, 
                         topic = c(10, 20, 30, 50, 80, 100),
                         frac_in_cell = 0.05){
  # prepare the right format of rownames
  if(!require(cisTopic)){
     devtools::install_github("aertslab/cisTopic")
  }
  library(cisTopic)
  rnames = data.table('region' = rownames(mtx))
  tmp = tidyr::separate(rnames, col = 'region', into = c('chr', 'start', 'end'))
  rnames = paste0(tmp$chr, ':', tmp$start, '-', tmp$end)
  rownames(mtx) = rnames
 
  ## reduce # of features to speed up
  mtx = 1 * (mtx > 0)
  rr = Matrix::rowMeans(mtx)
  mtx = mtx[rr >= frac_in_cell, ]
     
  cisTopicObject <- createcisTopicObject(mtx, project.name='scATAC')
  cisTopicObject <- runModels(cisTopicObject, topic = topic, seed = 987, nCores = nCores, 
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
    
    if(only.pos) res0 = res0[mean1 > mean2]
    
    res0 = res0[order(pv_adjust), ]
    res0 = res0[pv_adjust <= fdr]
    
    if(nrow(res0) > topn) res0 = res0[1:topn, ]
    res = rbind(res, res0)
  }
  return(res)
}
do_DA = cmpfun(do_DA)

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
  
  #return(enrich_list@result[enrich_list@result$qvalue <= qCutoff, ])
  return(enrich_list)
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
                       assay = assay, seed.use = 10, features = vFeatures)
  
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


normalize_gene_activities.corrected <- function (activity_matrices, cell_num_genes) 
{
  if (!is.list(activity_matrices)) {
    scores <- activity_matrices
    normalization_df <- data.frame(cell = colnames(activity_matrices), 
                                   cell_group = 1)
  }else {
    scores <- do.call(cbind, activity_matrices)
    normalization_df <- do.call(rbind, lapply(seq_along(activity_matrices), 
                                              function(x) {
                                                data.frame(cell = colnames(activity_matrices[[x]]), 
                                                           cell_group = rep(x, ncol(activity_matrices[[x]])))
                                              }))
  }
  scores <- scores[Matrix::rowSums(scores) != 0, Matrix::colSums(scores) != 
                     0]
  
  ## correct by adding following lines
  cell_num_genes = cell_num_genes[normalization_df$cell %in% colnames(scores)]
  normalization_df = subset(normalization_df, cell %in% colnames(scores))
  ##
  
  normalization_df$cell_group <- factor(normalization_df$cell_group)
  normalization_df$total_activity <- Matrix::colSums(scores)
  normalization_df$total_sites <- cell_num_genes[as.character(normalization_df$cell)]
  if (!is.list(activity_matrices)) {
    activity_model <- stats::lm(log(total_activity) ~ log(total_sites), 
                                data = normalization_df)
  } else {
    activity_model <- stats::lm(log(total_activity) ~ log(total_sites) * 
                                  cell_group, data = normalization_df)
  }
  normalization_df$fitted_curve <- exp(as.vector(predict(activity_model, 
                                                         type = "response")))
  size_factors <- log(normalization_df$fitted_curve)/mean(log(normalization_df$fitted_curve))
  size_factors <- Matrix::Diagonal(x = 1/size_factors)
  row.names(size_factors) <- normalization_df$cell
  colnames(size_factors) <- row.names(size_factors)
  
  scores <- Matrix::t(size_factors %*% Matrix::t(scores))
  scores@x <- pmin(1e+09, exp(scores@x) - 1)
  sum_activity_scores <- Matrix::colSums(scores)
  scale_factors <- Matrix::Diagonal(x = 1/sum_activity_scores)
  row.names(scale_factors) <- normalization_df$cell
  colnames(scale_factors) <- row.names(scale_factors)
  scores <- Matrix::t(scale_factors %*% Matrix::t(scores))
  
  if (!is.list(activity_matrices)) {
    rn = row.names(activity_matrices)
    cn = colnames(activity_matrices)
    rn = rn[rn %in% row.names(scores)]
    cn = cn[cn %in% colnames(scores)]
    
    ret <- scores[rn, cn]
  }
  else {
    ret <- lapply(activity_matrices, function(x) {
      rn = row.names(x)
      cn = colnames(x)
      rn = rn[rn %in% row.names(x)]
      cn = cn[cn %in% colnames(x)]
      scores[rn, cn]
    })
  }
  return(ret)
}
normalize_gene_activities.corrected = cmpfun(normalize_gene_activities.corrected)

# do cicero given a Seurat object, output gene activity score
doCicero_gascore <- function(seurat.obj, reduction = 'tsne', chr_sizes,
                            gene_ann, npc = 30, coaccess_thr = 0.25){
  ## gene_ann: the first four columns: chr, start, end, gene name
  set.seed(2019)
  library(cicero)
  mtx = GetAssayData(seurat.obj, slot = 'counts')
  # change rownames using _ to delimited
  rnames = rownames(mtx)
  new.rnames = sapply(rnames, function(x) unlist(strsplit(x, ','))[1])
  
  new.rnames = sapply(new.rnames, function(x) gsub('-', '_', x))
  rownames(mtx) <- new.rnames
  
  
  #dt = reshape2::melt(as.matrix(mtx), value.name = 'count')
  #dt = dt[dt$count > 0, ]
  dt = mefa4::Melt(mtx)
  rm(mtx)
  names(dt) = c('Peak', 'Cell', 'Count') 
  dt$Cell = as.character(dt$Cell)
  dt$Peak = as.character(dt$Peak)
  
  input_cds <- make_atac_cds(dt, binarize = T)
  rm(dt)
  input_cds <- detectGenes(input_cds)
  input_cds <- estimateSizeFactors(input_cds)
  
  if(reduction == 'tsne') {
    if(is.null(seurat.obj@reductions$tsne))
      seurat.obj <- RunTSNE(seurat.obj, dims = 1:npc, check_duplicates = F)
    redu.coords = seurat.obj@reductions$tsne@cell.embeddings
  }
  if(reduction == 'umap') {
    if(is.null(seurat.obj@reductions$umap))
      seurat.object <- RunUMAP(seurat.object, dims = 1:npc)
    redu.coords = seurat.obj@reductions$umap@cell.embeddings
  }
  
  #make the cell id consistet
  
  cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = redu.coords)
  
  ## get connections
  
  conns <- run_cicero(cicero_cds, chr_sizes)
  
  ## get cicero gene activity score
  names(gene_ann)[4] <- "gene"
  
  input_cds <- annotate_cds_by_site(input_cds, gene_ann)
  
  # generate unnormalized gene activity matrix
  unnorm_ga <- build_gene_activity_matrix(input_cds, conns)
  
  # make a list of num_genes_expressed
  num_genes <- pData(input_cds)$num_genes_expressed
  names(num_genes) <- row.names(pData(input_cds))
  
  # normalize
  cicero_gene_activities <- normalize_gene_activities.corrected(unnorm_ga, num_genes)
  
  # if you had two datasets to normalize, you would pass both:
  # num_genes should then include all cells from both sets
  #unnorm_ga2 <- unnorm_ga
  #cicero_gene_activities <- normalize_gene_activities(list(unnorm_ga, unnorm_ga2), num_genes)
  conns = data.table(conns)
  conns = conns[coaccess > coaccess_thr, ]
  res = list('conns' = conns, 'ga_score' = cicero_gene_activities)
  return(res)  
}
doCicero_gascore = cmpfun(doCicero_gascore)

# do cicero given a Seurat object, just return the connection 
doCicero_conn <- function(seurat.obj, reduction = 'tsne', 
                          chr_sizes, npc = 30, coaccess_thr = 0.25){
  ## gene_ann: the first four columns: chr, start, end, gene name
  set.seed(2019)
  library(cicero)
  mtx = GetAssayData(seurat.obj, slot = 'counts')
  # change rownames using _ to delimited
  rnames = rownames(mtx)
  new.rnames = sapply(rnames, function(x) unlist(strsplit(x, ','))[1])
  new.rnames = sapply(new.rnames, function(x) gsub('-', '_', x))
  rownames(mtx) <- new.rnames
  
  #dt = reshape2::melt(as.matrix(mtx), value.name = 'count')
  #dt = dt[dt$count > 0, ]
  dt = mefa4::Melt(mtx)
  rm(mtx)
  names(dt) = c('Peak', 'Cell', 'Count')
  dt$Cell = as.character(dt$Cell)
  dt$Peak = as.character(dt$Peak)
  input_cds <- make_atac_cds(dt, binarize = T)
  rm(dt)
  input_cds <- detectGenes(input_cds)
  input_cds <- estimateSizeFactors(input_cds)
  
  if(reduction == 'tsne') {
    if(is.null(seurat.obj@reductions$tsne))
      seurat.obj <- RunTSNE(seurat.obj, dims = 1:npc, check_duplicates = F)
    redu.coords = seurat.obj@reductions$tsne@cell.embeddings
  }
  if(reduction == 'umap') {
    if(is.null(seurat.obj@reductions$umap))
      seurat.object <- RunUMAP(seurat.object, dims = 1:npc)
    redu.coords = seurat.obj@reductions$umap@cell.embeddings
  }
  
  #make the cell id consistet
  
  cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = redu.coords)
  
  ## get connections
  
  conns <- run_cicero(cicero_cds, chr_sizes)
  conns = data.table(conns)
  conns = conns[coaccess > coaccess_thr, ]
  return(conns)
}
doCicero_conn = cmpfun(doCicero_conn)

## change rowname of zscores (tf name) from chromvar to be readable
readable_tf <- function(sele.zscores, GENOME_NAME){
  if(grepl(GENOME_NAME, pattern = 'hg', ignore.case = T)){
    rnames = rownames(sele.zscores)
    nnames = sapply(rnames, function(x) unlist(strsplit(x, '_'))[3])
    nnames1 = sapply(rnames, function(x) unlist(strsplit(x, '_'))[1])
    rownames(sele.zscores) = ifelse(grepl(nnames, pattern = 'LINE'), nnames1, nnames)
}else{
    rnames = rownames(sele.zscores)
    nnames = sapply(rnames, function(x) unlist(strsplit(x, '_'))[3])
    rownames(sele.zscores) = nnames
    sele.zscores = sele.zscores[!grepl(nnames, pattern = '^LINE'), ]
}
 return(sele.zscores)
}

## plot enriched tf for chromvar
plot_enrich_tf <- function(sele.zscores, bc_clusters,
                           up.qt = 0.95, low.qt = 0.05,
                           ndownsample = 1000){
    
    bc_clusters = data.table(bc_clusters)
    
    #downsample 
    set.seed(2019)
    bc_clusters.down = bc_clusters
    
    if(!is.null(ndownsample) & ndownsample < nrow(bc_clusters.down)) 
      bc_clusters.down = bc_clusters.down[sort(sample((1:nrow(bc_clusters.down)), ndownsample)), ]
    
    bc_clusters = bc_clusters.down
    rr = bc_clusters$barcode[bc_clusters$barcode %in% colnames(sele.zscores)]
    sele.zscores = sele.zscores[, rr]
    
    ann_column = data.frame('cluster' = as.integer(bc_clusters$cluster), 
                            'barcode' = bc_clusters$barcode,
                            stringsAsFactors = F)
    rownames(ann_column) = bc_clusters$barcode
    
    up_cut = quantile(sele.zscores, up.qt, na.rm = T)
    low_cut = quantile(sele.zscores, low.qt, na.rm = T)
    sele.zscores[is.na(sele.zscores)] = 0
    low_cut = min(0, low_cut)
    sele.zscores[sele.zscores > up_cut] = up_cut
    sele.zscores[sele.zscores < low_cut] = low_cut
    
    nc = length(unique(bc_clusters$cluster))
    getPalette = colorRampPalette(brewer.pal(9, "Paired"))
    if(nc >= 3) color_cluster = getPalette(nc)
    if(nc < 3) color_cluster = c("#A6CEE3", "#1F78B4", "#B2DF8A")[1:nc]
    names(color_cluster) = sort(unique(bc_clusters$cluster))
    
    ann_column = ann_column[order(ann_column$cluster), ]
    
    ann_column$barcode = NULL
   
    ann_colors = list('cluster' = color_cluster)
    
    ann_column$cluster = factor(ann_column$cluster,
                                levels = sort(unique(ann_column$cluster)))
    
    ph <- pheatmap::pheatmap(sele.zscores[, rownames(ann_column)], 
                             cluster_cols = F, cluster_rows = F, 
                             show_colnames = F, fontsize = 9,
                             annotation_col = ann_column, 
                             color = viridis(100),
                             annotation_colors = ann_colors, 
                             fontsize_row = 9)
    return(ph)
}


read_conf_file <- function(configure_user_file){
  
  system(paste0('grep = ', configure_user_file, " | grep -v ^# | awk -F= '{print $1}' | awk '{$1=$1;print}' >  /tmp/vrs.txt"))
  
  system(paste0('grep = ', configure_user_file, " | grep -v ^# | awk -F= '{print $2}' | awk -F# '{print $1}' | awk '{$1=$1;print}' >  /tmp/vls.txt"))
  
  vrs = readLines('/tmp/vrs.txt')
  vls = readLines('/tmp/vls.txt')
  for(i in 1:length(vrs)){
    assign(vrs[i], vls[i], envir = .GlobalEnv)
  }
}

## remove features active in less than min_frac_per_cluser in first class
## do downsample to max_cellPer_cluster
## support wilcox and t test only
## test-model: one-rest; control
runDiffMotifEnrich <- function(mtx_score, clusters, test = 'wilcox',
                        fdr = 0.01, topn = 10,
                        min_frac_per_cluster = 0.1,
                        max_cell_per_clust = 500,
                        test.mode = 'one-rest', control = NULL){
  set.seed(2020)
  clusters$cluster = as.character(clusters$cluster)
  cls = unique(clusters$cluster)
  res = NULL
  features = rownames(mtx_score)
  for(cluster0 in cls){
    bc0 = clusters[cluster == cluster0]$barcode
    mtx1 = mtx_score[, colnames(mtx_score) %in% bc0]
    if(test.mode == 'one-rest') {
      mtx2 = mtx_score[, !colnames(mtx_score) %in% bc0]
    }
    if(test.mode == 'control'){
      if(is.null(control)) stop('Should specific control cluster!')
      bc2 = clusters[cluster %in% control]$barcode
      mtx2 = mtx_score[, colnames(mtx_score) %in% bc2]
      if(cluster0 %in% control) next
    }
    mu1 = sapply(1:length(features), function(x) mean(mtx1[x, ]))
    mu2 = sapply(1:length(features), function(x) mean(mtx2[x, ]))
    
    pvs = rep(0.5, length(features))
    
    for(x in 1:length(features)){
      a1 = mtx1[x, ]
      a2 = mtx2[x, ]
      if(mu1[x] <= 0) next
      if(length(which(!is.na(a1))) < 2 || length(which(!is.na(a2))) < 2) next
      if(mean(a1 > 0) < min_frac_per_cluster) next
      if(!is.null(max_cell_per_clust)){
        if(length(a1) > max_cell_per_clust) a1 = sample(a1, max_cell_per_clust)
        if(length(a2) > max_cell_per_clust) a2 = sample(a2, max_cell_per_clust)
      }
      
      if(test == 'wilcox') pvs[x] = wilcox.test(a1, a2, alternative = 'greater')$p.value
      if(test == 't') pvs[x] = t.test(a1, a2, alternative = 'greater')$p.value
      
    }
    pvs.adj = p.adjust(pvs, method = 'fdr')
    res0 = data.table('feature' = features, 'cluster1' = cluster0,
                      'mean1' = mu1, 'mean0' = mu2, 
                      'pv' = pvs, 'pv_adjust' = pvs.adj)
    
    res0 = res0[mean1 > 0]
    res0 = res0[order(pv_adjust), ]
    res0 = res0[pv_adjust <= fdr]
    
    if(nrow(res0) > topn) res0 = res0[1:topn, ]
    res = rbind(res, res0)
  }
  return(res)
}
runDiffMotifEnrich = cmpfun(runDiffMotifEnrich)

## seurat_list: a list of seurat obj -- used since v1.5.1
run_integrateSeuObj <- function(seurat_list, integrate_by = 'VFACS',
                            top_variable_features = 5000, 
                            norm_by = 'tf-idf', nREDUCTION = 30,
                            reg.var = 'nFeature_ATAC',
                            resolution = 0.6, verbose = F){
  
  ## pool/integrate data into a seurat object
  if(integrate_by %in% c('seurat', 'cca', 'rpca')){
   
    # select features that are repeatedly variable across datasets for integration run PCA on each
    # dataset using these features
    features <- SelectIntegrationFeatures(object.list = seurat_list, 
                                          nfeatures = top_variable_features)
    seurat_list <- lapply(X = seurat_list, FUN = function(x) {
      x <- ScaleData(x, features = features, verbose = FALSE)
      x <- RunPCA(x, features = features, verbose = FALSE)
    })
    anchors <- FindIntegrationAnchors(object.list = seurat_list,
                                      anchor.features = features,
                                      reduction = ifelse(integrate_by == 'rpca', 'rpca', 'cca'))
    rm(seurat_list)
    seurat.merged <- IntegrateData(anchorset = anchors, dims = 1:nREDUCTION)
    DefaultAssay(seurat.merged) <- "integrated"
    seurat.merged <- ScaleData(seurat.merged, verbose = FALSE,
                            features = VariableFeatures(seurat.merged))
    seurat.merged <- RunPCA(seurat.merged, npcs = nREDUCTION, verbose = verbose)
  }else{
    message('Merge seurat objects...')
    seurat.merged = merge(seurat_list[[1]], seurat_list[-1])
  }
  
  if(integrate_by == 'pool') seurat.merged <- regress_on_pca(seurat.merged, 'sampleName')
  
  if(integrate_by == 'VFACS'){
    ## cluster and then reselect features
    ## variable features across clusters
    message('Workin on merged object ...')
    
    if(class(seurat.merged[['ATAC']]) == 'Assay5'){
      seurat.merged[['ATAC']]$data = TF_IDF(seurat.merged[['ATAC']]$counts)
    }else{
      seurat.merged[['ATAC']]@data = TF_IDF(seurat.merged[['ATAC']]@counts)
    }
    
    seurat.merged <- FindVariableFeatures(seurat.merged, nfeatures = top_variable_features)
    seurat.merged <- ScaleData(seurat.merged)
    seurat.merged <- RunPCA(seurat.merged, npcs = nREDUCTION, verbose = verbose)
    seurat.merged <- FindNeighbors(seurat.merged, dims = 1:nREDUCTION, reduction = 'pca', 
                                verbose = verbose)
    seurat.merged <- FindClusters(seurat.merged, resolution = resolution, verbose = verbose)
    clusters = as.character(seurat.merged$seurat_clusters)
    
    if(class(seurat.merged[['ATAC']]) == 'Assay5') {
      mtx = seurat.merged[['ATAC']]$counts
    }else{
      mtx = seurat.merged[['ATAC']]@counts
    }
    mtx_by_cls <- sapply(unique(clusters), function(x) {
      
      cl_data <- mtx[, clusters == x]
      
      Matrix::rowMeans(cl_data > 0)
      
    })
    mtx_by_cls.norm <- edgeR::cpm(mtx_by_cls, log = T, prior.count = 1)
    sds = sapply(1:nrow(mtx_by_cls.norm), function(x) sd(mtx_by_cls.norm[x, ]))
    names(sds) = rownames(mtx_by_cls.norm)
    sele.features = names(which(sds >= sort(sds, decreasing = T)[top_variable_features]))
    mtx0 = mtx[sele.features, ]
    mtx0.norm = TF_IDF(mtx0)
    
    tmp <- mtx0[setdiff(rownames(mtx0), sele.features), ]
    data0 <- rbind(mtx0.norm, tmp)
    seurat.merged[['ATAC']]@data = data0[rownames(mtx0), ]
    
    VariableFeatures(seurat.merged) <- sele.features
    seurat.merged <- RunPCA(seurat.merged, npcs = nREDUCTION, verbose = verbose)
    seurat.merged <- regress_on_pca(seurat.merged, reg.var = reg.var)
    if(class(seurat.merged[['ATAC']]) == 'Assay5') {
      seurat.merged[['ATAC']]$data = TF_IDF(mtx)
    }else{
      seurat.merged[['ATAC']]@data <- TF_IDF(mtx)
    }
    
  }
  
  if(integrate_by %in% c('rlsi', 'signac')){
    library(Signac)
    ## process each sample
    for(sample0 in names(seurat_list)){
      seurat_list[[sample0]] = FindTopFeatures(seurat_list[[sample0]],
                                               min.cutoff = as.integer(0.01 * ncol(seurat_list[[sample0]])))
      seurat_list[[sample0]] = RunTFIDF(seurat_list[[sample0]])
      seurat_list[[sample0]] = RunSVD(seurat_list[[sample0]])
      
    }
    ## merge 
    seurat.merged = merge(seurat_list[[1]], seurat_list[-1])
    seurat.merged = FindTopFeatures(seurat.merged, min.cutoff = 100)
    
    ## select anchor features
    if(class(seurat.merged[['ATAC']]) == 'Assay5'){
      mtx = seurat.merged[['ATAC']]$counts
    }else{
      mtx = seurat.merged[['ATAC']]@counts
    }
  
    sele_features = rownames(mtx)
    sIDs = seurat.merged$sampleName
    mtx_pbulk <- sapply(unique(sIDs), function(x){
      rowMeans(mtx[, sIDs == x] >0)
    })
    rmaxs = apply(mtx_pbulk, 1, max)
    filtered_peaks1 <- names(which(rmaxs < 0.02)) # filter variable features
    sele_features = setdiff(sele_features, filtered_peaks1)
    sele_features = intersect(sele_features, VariableFeatures(seurat.merged))
    
    VariableFeatures(seurat.merged) <- sele_features
    seurat.merged = RunTFIDF(seurat.merged)
    seurat.merged = RunSVD(seurat.merged)
    seurat.merged = RunUMAP(seurat.merged, reduction = 'lsi', dims = 2:nREDUCTION)
    
    ## integration 
    integration.anchors <- FindIntegrationAnchors(
      object.list = seurat_list,
      anchor.features = sele_features,
      reduction = "rlsi",
      dims = 2:nREDUCTION
    )
    
    # integrate LSI embeddings
    seurat.merged <- IntegrateEmbeddings(
      anchorset = integration.anchors,
      reductions = seurat.merged[["lsi"]],
      new.reduction.name = "integrated_lsi",
      dims.to.integrate = 2:nREDUCTION,
      k.weight = 30
    )
    
    seurat.merged <- RunUMAP(seurat.merged, reduction = "integrated_lsi", dims = 2:nREDUCTION)
    seurat.merged = FindNeighbors(seurat.merged, reduction = 'integrated_lsi', 
                                  dims = 2:nREDUCTION, verbose = verbose)
    
  }
  
  if(integrate_by %in% c('seurat', 'cca', 'rpca',
                         'pool', 'VFACS')) {
    seurat.merged <- RunUMAP(seurat.merged, reduction = "pca", 
                        verbose = verbose, dims = 1:nREDUCTION)
    seurat.merged = FindNeighbors(seurat.merged, reduction = 'pca', 
                                  dims = 1:nREDUCTION, verbose = verbose)
    
  }
  
  if(integrate_by == 'harmony'){
    
    seurat.merged <- harmony::RunHarmony(seurat.merged, c("sampleName"), 
                                         assay.use = 'ATAC')
    
    seurat.merged <- seurat.merged %>% 
      RunUMAP(reduction = "harmony", dims = 1:nREDUCTION, verbose = verbose)
    
    seurat.merged = FindNeighbors(seurat.merged, reduction = 'harmony', 
                                  dims = 1:nREDUCTION, verbose = verbose)
    
  }
  
  ## clustering on the integrated data
  seurat.merged = FindClusters(seurat.merged, resolution = resolution, verbose = verbose)
  seurat.merged$active_clusters = seurat.merged$seurat_clusters
  
  return(seurat.merged)
}
run_integrateSeuObj = cmpfun(run_integrateSeuObj)   


## mtx_list: a list of matrix  -- older version
run_integration <- function(mtx_list, integrate_by = 'VFACS',
                            top_variable_features = 5000, 
                            norm_by = 'tf-idf', nREDUCTION = 30,
                            minFrac_in_cell = 0.01, min_depth = 1000,
                            max_depth = 50000, reg.var = 'nCount_ATAC',
                            anchor.features = 2000,
                            resolution = 0.6, verbose = F){
  nsample = length(mtx_list)
  sampleNames = names(mtx_list)
  if(is.null(sampleNames)) sampleNames = paste0('sample', 1:nsample)
  names(mtx_list) = sampleNames
  seu.all <- list()
  for(sample0 in sampleNames){
    # filter each mtx
    mtx_list[[sample0]] <- filterMat(mtx_list[[sample0]], minFrac_in_cell = minFrac_in_cell,
                                     min_depth = min_depth, max_depth = max_depth)
    if(integrate_by %in% c('seurat')){
      # create a seurat obj for each sample
      seurat.obj = runSeurat_Atac(mtx_list[[sample0]], npc = nREDUCTION, norm_by = norm_by, 
                                  top_variable_features = top_variable_features, 
                                  reg.var = reg.var)
      
      seurat.obj$sample = sample0
      
      seu.all[[sample0]] = seurat.obj
    }
  }
  
  ## pool/integrate data into a seurat object
  if(integrate_by == 'seurat'){
    
    seurat.obj <- FindIntegrationAnchors(object.list = seu.all, 
                                         anchor.features = anchor.features)
    rm(seu.all)
    seurat.obj <- IntegrateData(anchorset = seurat.obj, dims = 1:nREDUCTION)
    DefaultAssay(seurat.obj) <- "integrated"
    seurat.obj <- ScaleData(seurat.obj, verbose = FALSE,
                            features = VariableFeatures(seurat.obj))
    seurat.obj <- RunPCA(seurat.obj, npcs = nREDUCTION, verbose = verbose)
  }else{
    # pool the matrxi first with union features
    nf = sapply(mtx_list, nrow)
    nc = sapply(mtx_list, ncol)
    
    umtx <- cBind_union_features(mtx_list)
    
    rm(mtx_list)
    seurat.obj = runSeurat_Atac(umtx, npc = nREDUCTION, norm_by = norm_by, 
                                top_variable_features = top_variable_features, 
                                reg.var = reg.var)
    seurat.obj$sample = rep(sampleNames, nc)
    
  }
  
  if(integrate_by == 'pool') seurat.obj <- regress_on_pca(seurat.obj, 'sample')
  
  if(integrate_by == 'VFACS'){
    ## cluster and then reselect features
    ## variable features across clusters
    seurat.obj <- FindNeighbors(seurat.obj, dims = 1:nREDUCTION, reduction = 'pca', 
                                verbose = verbose)
    seurat.obj <- FindClusters(seurat.obj, resolution = resolution, verbose = verbose)
    clusters = as.character(seurat.obj$seurat_clusters)
    
    if(class(seurat.obj[['ATAC']]) == 'Assay5'){
      mtx = seurat.obj[['ATAC']]$counts
    }else{
      mtx = seurat.obj[['ATAC']]@counts
    }
      
    mtx_by_cls <- sapply(unique(clusters), function(x) {
      
      cl_data <- mtx[, clusters == x]
      
      Matrix::rowMeans(cl_data > 0)
      
    })
    mtx_by_cls.norm <- edgeR::cpm(mtx_by_cls, log = T, prior.count = 1)
    sds = sapply(1:nrow(mtx_by_cls.norm), function(x) sd(mtx_by_cls.norm[x, ]))
    names(sds) = rownames(mtx_by_cls.norm)
    sele.features = names(which(sds >= sort(sds, decreasing = T)[top_variable_features]))
    mtx0 = mtx[sele.features, ]
    mtx0.norm = TF_IDF(mtx0)
    
    if(class(seurat.obj[['ATAC']]) == 'Assay5') {
      seurat.obj[['ATAC']]$data[sele.features, ] <- mtx0.norm
    }else{
      seurat.obj[['ATAC']]@data[sele.features, ] <- mtx0.norm
    }
    
    VariableFeatures(seurat.obj) <- sele.features
    seurat.obj <- RunPCA(seurat.obj, dims = 1:nREDUCTION, verbose = verbose)
    seurat.obj <- regress_on_pca(seurat.obj, reg.var = reg.var)
    seurat.obj <- FindNeighbors(seurat.obj, verbose = verbose, 
                                dims = 1:nREDUCTION, reduction = 'pca')
    seurat.obj <- FindClusters(seurat.obj, verbose = verbose, resolution = resolution)
    
  }
  
  seurat.obj <- RunUMAP(seurat.obj, reduction = "pca", 
                        verbose = verbose, dims = 1:nREDUCTION)
  
  
  if(integrate_by == 'harmony'){
    
    seurat.obj <- harmony::RunHarmony(seurat.obj, c("sample"), assay.use = 'ATAC')
    
    seurat.obj <- seurat.obj %>% 
      RunUMAP(reduction = "harmony", dims = 1:nREDUCTION, verbose = verbose) 
  }
  
  ## clustering on the integrated data
  ## seurat implemented louvain algorithm
  redm = ifelse(integrate_by == 'harmony', 'harmony', 'pca')
  seurat.obj = FindNeighbors(seurat.obj, reduction = redm, 
                             dims = 1:nREDUCTION, verbose = verbose)
  
  seurat.obj = FindClusters(seurat.obj, resolution = resolution, verbose = verbose)
  seurat.obj$active_clusters = seurat.obj$seurat_clusters
  
  
  return(seurat.obj)
}
run_integration = cmpfun(run_integration)   



# Find doublets
FindDoublets_Atac <- function(seurat.atac, PCs = 1:50, 
                         exp_rate = 0.02, sct = FALSE){
  # sct--do SCTransform or not
  require('DoubletFinder') 
  if(!sct){
    if(class(seurat.atac[['ATAC']]) == 'Assay5') {
      seurat.rna = CreateSeuratObject(seurat.atac[['ATAC']]$counts)
    }else{
      seurat.rna = CreateSeuratObject(seurat.atac[['ATAC']]@counts)
    }
    seurat.rna = NormalizeData(seurat.rna)
    seurat.rna = FindVariableFeatures(seurat.rna)
    VariableFeatures(seurat.rna) <- VariableFeatures(seurat.atac)
    seurat.rna = ScaleData(seurat.rna)
    seurat.rna = RunPCA(seurat.rna, npcs = max(PCs), verbose = F)
    seurat.rna@reductions$pca@cell.embeddings <- seurat.atac@reductions$pca@cell.embeddings
    seurat.rna@reductions$pca@feature.loadings <- seurat.atac@reductions$pca@feature.loadings
    seurat.rna$seurat_clusters = seurat.atac$seurat_clusters
  }else{
    seurat.rna = seurat.atac
    seurat.rna[['RNA']] <- seurat.atac[['ATAC']]
  }
  
  ## pK identification
  sweep.res.list <- paramSweep(seurat.rna, PCs = PCs, sct = sct)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  ## Homotypic Doublet proportion Estimate
  annotations <- seurat.rna@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(exp_rate * length(seurat.rna$seurat_clusters))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies
  seurat.rna <- doubletFinder(seurat.rna, PCs = PCs, pN = 0.25,
                                 pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, 
                                 sct = sct)
  
  seurat.rna <- doubletFinder(seurat.rna, PCs = PCs, pN = 0.25, 
                                 pK = 0.09, nExp = nExp_poi.adj,
                                 reuse.pANN = paste0("pANN_0.25_0.09_", nExp_poi), 
                                 sct = sct)
  doublet_var = paste0('DF.classifications_0.25_0.09_', nExp_poi.adj)
  seurat.rna[['Doublet_Singlet']] = seurat.rna[[doublet_var]]
  
  mnames = names(seurat.rna@meta.data)
  seurat.rna@meta.data[, grep(mnames, pattern = '0.25_0.09')] <- NULL
  seurat.atac$Doublet_Singlet = seurat.rna$Doublet_Singlet
  return(seurat.atac)
}
FindDoublets_Atac = cmpfun(FindDoublets_Atac)

# generate gene activity matrix for labelTransfer
generate_gene_cisActivity <- function(gene_ann, mtx, include_body = T,
                                      dist_to_tss = 2000){
  ## generating gene cis activity score
  ## input: gene_ann (data table with chr, gene_start, gene_end, strand, gene_name), 
  ##        and atac matrix file
  
  ## 1. select gene up/down-stream regions (promoter with/without gene_body) ##
  
  
  if(!include_body){
    gene_ann[, 'start' := ifelse(strand == '+', gene_start - dist_to_tss, gene_end - dist_to_tss)]
    gene_ann[, 'end' := ifelse(strand == '+', gene_start + dist_to_tss, gene_end + dist_to_tss)]
  }else{
    gene_ann[, 'start' := ifelse(strand == '+', gene_start - dist_to_tss, gene_start)]
    gene_ann[, 'end' := ifelse(strand == '+', gene_end, gene_end + dist_to_tss)]
    
  }
  
  gene_ann = subset(gene_ann, select = c(chr, start, end, gene_name))
  
  
  ## 2. read mtx file ##
  
  rnames = rownames(mtx)
  chrs = sapply(rnames, function(x) unlist(strsplit(x, '-'))[1])
  starts = sapply(rnames, function(x) unlist(strsplit(x, '-'))[2])
  ends = sapply(rnames, function(x) unlist(strsplit(x, '-'))[3])
  
  peaks = data.table('chr' = chrs, 'start' = as.integer(starts), 
                     'end' = as.integer(ends))
  setkey(peaks, chr, start, end)
  peaks[, 'pname' := paste(chr, start, end, sep = '-')]
  over.ids = foverlaps(gene_ann, peaks, by.x = c('chr', 'start', 'end'),
                       by.y = c('chr', 'start', 'end'), which = T)
  over.ids[, 'gene_name' := gene_ann[xid, gene_name]]
  over.ids[, 'pname' := peaks[yid, pname]]
  over.ids = over.ids[complete.cases(over.ids)]
  
  
  smtx = sparseMatrix(i = over.ids$xid, j = over.ids$yid,
                      dimnames = list(gene_ann$gene_name[1:max(over.ids$xid)],
                                      peaks$pname[1:max(over.ids$yid)]))
  
  mtx = mtx[rownames(mtx) %in% colnames(smtx), ]
  smtx = smtx[, rownames(mtx)]
  
  activity.matrix = smtx %*% mtx
  rs = Matrix::rowSums(activity.matrix)
  activity.matrix = activity.matrix[rs > 10, ]
  
  return(activity.matrix)
  
}

generate_gene_cisActivity = cmpfun(generate_gene_cisActivity)


labelTransfer_R <- function(seurat.atac, seurat.rna, gene_ann,
                            rna_ann_var = 'Cell_Type', include_genebody = T){
  if(class(seurat.atac[['ATAC']]) == 'Assay5'){
    atac.mtx = seurat.atac[['ATAC']]$counts
  }else{
    atac.mtx = seurat.atac[['ATAC']]@counts
  }
  rn = rownames(atac.mtx)
  rownames(atac.mtx) <- sapply(rn, function(x) unlist(strsplit(x, ','))[1])
  activity.matrix = generate_gene_cisActivity(gene_ann = gene_ann,
                                              atac.mtx, 
                                              include_body = include_genebody)
  
  seurat.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
  
  DefaultAssay(seurat.atac) <- "ACTIVITY"
  seurat.atac <- NormalizeData(seurat.atac)
  seurat.atac <- ScaleData(seurat.atac, features = rownames(seurat.atac))
  seurat.atac <- FindVariableFeatures(seurat.atac)
  
  DefaultAssay(seurat.atac) <- "ATAC"
  seurat.atac$tech = 'ATAC'
  seurat.rna$tech = 'RNA'
  
  ## transfer label 
  genes4anchors = VariableFeatures(object = seurat.rna)
  #genes4anchors = NULL
  transfer.anchors <- FindTransferAnchors(reference = seurat.rna,
                                          query = seurat.atac,
                                          features = genes4anchors,
                                          reference.assay = "RNA",
                                          query.assay = "ACTIVITY",
                                          reduction = "cca",
                                          k.anchor = 5)
  
  
  celltype.predictions <- TransferData(anchorset = transfer.anchors,
                                       refdata = seurat.rna[[rna_ann_var]][, 1],
                                       weight.reduction = seurat.atac[["pca"]],
                                       dims = 1:ncol(seurat.atac[["pca"]]),
                                       k.weight = 50)
  celltype.predictions = subset(celltype.predictions, select = c('predicted.id', 'prediction.score.max'))
  names(celltype.predictions)[1] = 'Predicted_Cell_Type'
  seurat.atac <- AddMetaData(seurat.atac, metadata = celltype.predictions)
  rm(transfer.anchors)
  
  seurat.atac[["ACTIVITY"]] <- NULL ## don't save activity assay
  return(seurat.atac)
}
labelTransfer_R = cmpfun(labelTransfer_R)
