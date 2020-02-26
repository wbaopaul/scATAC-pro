source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local('dsAnalysis_utilities.R')

library(parallel)


args = commandArgs(T)
mtx_file = args[1]
genome_name = args[2]
output_dir = args[3]
norm_by = args[4]

mtx = read_mtx_scATACpro(mtx_file)
#mtx = filterMat(mtx)
genomeName = 'BSgenome.Hsapiens.UCSC.hg38'
if(grepl(genome_name, pattern = '38'))genomeName = 'BSgenome.Hsapiens.UCSC.hg38'
if(grepl(genome_name, pattern = '19') || grepl(genome_name, pattern = '37'))genomeName = 'BSgenome.Hsapiens.UCSC.hg19'
if(grepl(genome_name, pattern = 'mm9'))genomeName = 'BSgenome.Mmusculus.UCSC.mm9'
if(grepl(genome_name, pattern = 'mm10'))genomeName = 'BSgenome.Mmusculus.UCSC.mm10'

ncore = detectCores()

if(T){
## select variable features first
      seurat.obj = CreateSeuratObject(mtx, project = 'scATAC', assay = 'ATAC',
                                  names.delim = '-')

      if(norm_by == 'log') seurat.obj@assays$ATAC@data = log1p(seurat.obj@assays$ATAC@counts)
      if(norm_by == 'tf-idf') seurat.obj@assays$ATAC@data = TF.IDF(seurat.obj@assays$ATAC@counts)
      
      seurat.obj <- FindVariableFeatures(object = seurat.obj,
                                         selection.method = 'vst',
                                         nfeatures = floor(nrow(mtx) * 0.3))
      vFeatures = VariableFeatures(seurat.obj)
      rm(seurat.obj)
      rnames = rownames(mtx)
      mtx = mtx[rnames %in% vFeatures, ]
}


chromVar.obj = run_chromVAR(mtx, genomeName, max(1, ncore - 1))
saveRDS(chromVar.obj, file = paste0(output_dir, '/chromVar_obj.rds'))


## plot heatmap ####
library(BiocParallel)
register(SerialParam())

# Do DA/DE with one cluster vs the rest clusters
# clusters are the data frame with <barcode> <cluster>
do_DA_motif <- function(mtx_score, clusters, test = 'wilcox', 
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
    
    pvs = rep(0.5, length(features))
    
    for(x in 1:length(features)){
      a1 = mtx1[x, ]
      a2 = mtx2[x, ]
      if(length(which(!is.na(a1))) < 2 || length(which(!is.na(a2))) < 2) next
      pvs[x] = wilcox.test(a1, a2, alternative = 'greater')$p.value
    }
    
    pvs.adj = p.adjust(pvs, method = 'fdr')
    res0 = data.table('feature' = features, 'cluster' = cluster0,
                      'mean1' = mu1, 'mean2' = mu2,
                      'pv' = pvs, 'pv_adjust' = pvs.adj)
    
    
    res0 = res0[order(pv_adjust), ]
    res0 = res0[pv_adjust <= fdr]
    
    if(nrow(res0) > topn) res0 = res0[1:topn, ]
    res = rbind(res, res0)
  }
  return(res)
}

seurat_file <- paste0(output_dir, '/seurat_obj.rds')
if(file.exists(seurat_file)){
  ss = readRDS(seurat_file)
  metaData = ss@meta.data
  rm(ss)
  if(file.exists(paste0(output_dir, '/chromVar_obj.rds'))){
    #chromVar.obj = readRDS(paste0(output_dir, '/chromVar_obj.rds'))
    
    dev = deviations(chromVar.obj)
    da.res = do_DA_motif(dev, 
                   clusters = data.table('barcode' = rownames(metaData),
                                         'cluster' = metaData$active_clusters),
                   topn = 10)
    rm(dev)
    write.csv(da.res, file = paste0(output_dir, '/differential_TF_motif_enriched_in_clusters.txt'), 
              quote = F, row.names = F )
    
    
    ## plot enriched TFs in heatmap
    sele.tfs = da.res$feature
    #zscores = chromVar.obj@assays$data$z
    zscores = deviationScores(chromVar.obj)
    sele.zscores = zscores[sele.tfs, ]
    
    
    # change tf name to be more readable
    if(grepl(genome_name, pattern = 'hg', ignore.case = T)){
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
    metaData$active_clusters = as.character(metaData$active_clusters)
    metaData = data.table(metaData, keep.rownames = T)
    setkey(metaData, active_clusters)
    
    rr = metaData$rn[metaData$rn %in% colnames(sele.zscores)]
    sele.zscores = sele.zscores[, rr]
    
    
    sele.zscores = sele.zscores[!duplicated(sele.zscores), ]
    
    ann_col = data.frame('cluster' = metaData$active_clusters)
    rownames(ann_col) = metaData$rn
    
    up_cut = quantile(sele.zscores, 0.95, na.rm = T)
    low_cut = quantile(sele.zscores, 0.05, na.rm = T)
    sele.zscores[is.na(sele.zscores)] = 0
    low_cut = min(0, low_cut)
    sele.zscores[sele.zscores > up_cut] = up_cut
    sele.zscores[sele.zscores < low_cut] = low_cut
    
    cluster = brewer.pal(n=length(unique(metaData$active_clusters)), name = 'Paired')
    names(cluster) = sort(unique(metaData$active_clusters))
    ann_colors = list('cluster' = cluster)
    
    # resample to reduce memory used
    set.seed(2019)
    rids = sort(sample(1:ncol(sele.zscores), floor(ncol(sele.zscores)/6)))
    ann_col0 = data.frame(ann_col[rids, ])
    rownames(ann_col0) = colnames(sele.zscores)[rids]
    mtx0 = sele.zscores[, rids]
    names(ann_col0) = 'cluster'
    ph <- pheatmap::pheatmap(mtx0, cluster_cols = F,
                             cluster_rows = T, show_colnames = F, fontsize = 13,
                             annotation_col = ann_col0, color = viridis(100),
                             annotation_colors = ann_colors, fontsize_row = 9)
    
   
    pfname = paste0(params$output_dir, '/heatmap_motif_enrich.eps')
    #postscript(file = pfname, width = 9, height = 12)
    
    ggsave(ph, filename = pfname, device = 'eps', height = 12,
           width = 9)
    #dev.off()
    
    
  }
  
}



