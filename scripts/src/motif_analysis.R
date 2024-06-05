source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local('dsAnalysis_utilities.R')

library(parallel)



args = commandArgs(T)
mtx_file = args[1] ## could be the path of matrix.mtx or rds file for matrix or seurat obj
genome_name = args[2]
output_dir = args[3]
norm_by = args[4]

if(grepl(mtx_file, pattern = '.rds', fix = T)) {
    input.obj = readRDS(mtx_file)
    itype = class(input.obj)
    if(any(itype == 'Seurat')){
      if(class(input.obj[['ATAC']]) == 'Assay5'){
        mtx = input.obj[['ATAC']]$counts
      }else{
        mtx = input.obj[['ATAC']]@counts
      }
    }else{
        mtx = input.obj
        rm(input.obj)
    }
    rnames = rownames(mtx)
    new.rnames = sapply(rnames, function(x) unlist(strsplit(x, ','))[1])
    new.rnames = sapply(new.rnames, function(x) gsub('_', '-', x))
    names(new.rnames) = NULL
    rownames(mtx) <- new.rnames
}else{
    mtx = read_mtx_scATACpro(mtx_file)
}

#mtx = filterMat(mtx)
genomeName = 'BSgenome.Hsapiens.UCSC.hg38'
if(grepl(genome_name, pattern = '38'))genomeName = 'BSgenome.Hsapiens.UCSC.hg38'
if(grepl(genome_name, pattern = '19') || grepl(genome_name, pattern = '37'))genomeName = 'BSgenome.Hsapiens.UCSC.hg19'
if(grepl(genome_name, pattern = 'mm9'))genomeName = 'BSgenome.Mmusculus.UCSC.mm9'
if(grepl(genome_name, pattern = 'mm10'))genomeName = 'BSgenome.Mmusculus.UCSC.mm10'

ncore = detectCores()

if(all(itype != 'Seurat')) {
    seurat.obj = CreateSeuratObject(mtx, project = 'scATAC', assay = 'ATAC',
                          names.delim = '-')
    if(class(seurat.obj[['ATAC']]) == 'Assay5'){
      if(norm_by == 'log') seurat.obj[['ATAC']]$data = log1p(seurat.obj[['ATAC']]$counts)
      if(norm_by == 'tf-idf') seurat.obj[['ATAC']]$data = TF_IDF(seurat.obj[['ATAC']]$counts)
    }else{
      if(norm_by == 'log') seurat.obj[['ATAC']]@data = log1p(seurat.obj[['ATAC']]@counts)
      if(norm_by == 'tf-idf') seurat.obj[['ATAC']]@data = TF_IDF(seurat.obj[['ATAC']]@counts)
    }
    
    seurat.obj <- FindVariableFeatures(object = seurat.obj,
                                 selection.method = 'vst',
                                 nfeatures = floor(nrow(mtx) * 0.4))
    vFeatures = VariableFeatures(seurat.obj)
    ## further filter peaks
    rs = Matrix::rowSums(mtx > 0)
    filter.pks = names(which(rs > (0.005 * ncol(seurat.obj))))
    vFeatures = intersect(vFeatures, filter.pks)

    rm(seurat.obj)
    rnames = rownames(mtx)
    mtx = mtx[rnames %in% vFeatures, ]
}else{
    seurat.obj = input.obj
    rm(input.obj)
}



chromVar.obj = run_chromVAR(mtx, genomeName, max(1, ncore - 1))
seurat_file <- paste0(output_dir, '/seurat_obj.rds') ## default seurat obj for the project, not used if the input is an seurat obj
if(any(itype == 'Seurat')) {
    output_dir = dirname(mtx_file)
}
saveRDS(chromVar.obj, file = paste0(output_dir, '/chromVar_obj.rds'))

## plot enrichment heatmap ####
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

if(file.exists(seurat_file) | any(itype == 'Seurat')){
    if(all(itype != 'Seurat')){
        seurat.obj = readRDS(seurat_file)
    }

    metaData = seurat.obj@meta.data
    rm(seurat.obj)

    dev = deviations(chromVar.obj)
    da.res = do_DA_motif(dev, 
                   clusters = data.table('barcode' = rownames(metaData),
                                         'cluster' = metaData$active_clusters),
                   topn = 10)
    rm(dev)
    write.table(da.res, file = paste0(output_dir, '/differential_TF_motif_enriched_in_clusters.tsv'), 
              quote = F, sep = '\t', row.names = F )


    ## plot enriched TFs in heatmap
    sele.tfs = da.res$feature
    #zscores = chromVar.obj@assays$data$z
    zscores = deviationScores(chromVar.obj)
    sele.zscores = zscores[sele.tfs, ]

    # change tf name to be more readable
    sele.zscores = readable_tf(sele.zscores, genome_name)    

    metaData$active_clusters = as.character(metaData$active_clusters)

    bc_clusters = data.table('barcode' = rownames(metaData),
                             'cluster' = metaData$active_clusters)  

    ph <- plot_enrich_tf(sele.zscores, bc_clusters) 
    pfname = paste0(output_dir, '/heatmap_motif_enrich.eps')

    ggsave(ph, filename = pfname, device = 'eps', height = 12,
           width = 9)

}



