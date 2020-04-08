source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local('dsAnalysis_utilities.R')

library(ggplot2)
library(parallel)
args = commandArgs(T)
mtx_files = args[1]
k = (args[2])
output_dir = args[3]
genome_name = args[4]
tss_path = args[5]
norm_by = args[6]
REDUCTION = args[7]
nREDUCTION = as.integer(args[8])
top_variable_features = as.numeric(args[9])
integrate_by = args[10]

message(paste('Integrate by', integrate_by))

mtx_files = unlist(strsplit(mtx_files, ','))

tss_ann <- fread(tss_path, header = F)
#names(tss_ann)[c(1:4,7)] <- c('chr', 'start', 'end', 'gene_name', 'gene_type')
#tss_ann <- tss_ann[gene_type %in% c('miRNA', 'lincRNA', 'protein_coding'), ]
names(tss_ann)[c(1:4)] <- c('chr', 'start', 'end', 'gene_name')


## do seurat individually
seu.all = mtx.all = list()
len = length(mtx_files)
for(i in 1:length(mtx_files)){
    file0 = mtx_files[i]
    mtx = read_mtx_scATACpro(file0)
    rs = Matrix::rowMeans(mtx > 0)
    mtx = mtx[rs > 0.01, ]
    mtx = assignGene2Peak(mtx, tss_ann)
    #colnames(mtx) = paste0('sample', i, '_', colnames(mtx))
    if(integrate_by != 'seurat') mtx.all[[i]] = mtx
    nveg0 = ifelse(top_variable_features > 1, top_variable_features, floor(top_variable_features)*nrow(mtx))
    nveg = ifelse(nveg0 < nrow(mtx)/2, nveg0, floor(nrow(mtx)/2))
    seurat.obj = doBasicSeurat_new(mtx, npc = nREDUCTION, norm_by = norm_by, 
                                    top_variable_features = nveg, 
                                       reg.var = 'nCount_ATAC')
    
    seurat.obj$sample = paste0('sample', i)
   
   seu.all[[i]] = seurat.obj
   rm(seurat.obj, mtx)
}

if(integrate_by == 'seurat'){
    seurat.obj <- FindIntegrationAnchors(object.list = seu.all)
    rm(seu.all)
    seurat.obj <- IntegrateData(anchorset = seurat.obj, dims = 1:nREDUCTION)
    DefaultAssay(seurat.obj) <- "integrated"
    seurat.obj <- ScaleData(seurat.obj, verbose = FALSE,
                         features = VariableFeatures(seurat.obj))
    seurat.obj <- RunPCA(seurat.obj, npcs = nREDUCTION, verbose = FALSE)
}else{
    ## use pool and regress method
    nf = sapply(mtx.all, nrow)
    nc = sapply(mtx.all, ncol)
    if(length(unique(nf)) == 1) {
      umtx = do.cbind(mtx.all)  
    }else{
      umtx <- cBind_union_features(mtx.all)
    }
    rm(mtx.all)
    nveg0 = ifelse(top_variable_features > 1, top_variable_features, floor(top_variable_features)*nrow(umtx))
    nveg = ifelse(nveg0 < nrow(umtx)/2, nveg0, floor(nrow(umtx)/2))
    seurat.obj = doBasicSeurat_new(umtx, npc = nREDUCTION, norm_by = norm_by, 
                                    top_variable_features = nveg, 
                                       reg.var = 'nCount_ATAC')
    seurat.obj$sample = paste0('sample', rep(c(1:length(nc)), nc))
}

if(integrate_by == 'pool') seurat.obj <- regress_on_pca(seurat.obj, 'sample')

if(integrate_by == 'VFACS'){
    ## cluster and then reselect features
    ## variable features across clusters
    seurat.obj <- FindNeighbors(seurat.obj, dims = 1:nREDUCTION, reduction = 'pca')
    seurat.obj <- FindClusters(seurat.obj, resl = 0.6)
    clusters = as.character(seurat.obj$seurat_clusters)
    mtx = seurat.obj@assays$ATAC@counts
      mtx_by_cls <- sapply(unique(clusters), function(x) {
        
        cl_data <- mtx[, clusters == x]
        
        Matrix::rowMeans(cl_data)
        
      })
      mtx_by_cls.norm <- edgeR::cpm(mtx_by_cls, log = T, prior.count = 1)
      sds = sapply(1:nrow(mtx_by_cls.norm), function(x) sd(mtx_by_cls.norm[x, ]))
      names(sds) = rownames(mtx_by_cls.norm)
      sele.features = names(which(sds >= sort(sds, decreasing = T)[nveg]))
      mtx0 = mtx[sele.features, ]
      mtx0.norm = Seurat::TF.IDF(mtx0)
      seurat.obj@assays$ATAC@data[sele.features, ] <- mtx0.norm
      VariableFeatures(seurat.obj) <- sele.features
      seurat.obj <- RunPCA(seurat.obj, dims = 1:nReduction, verbose = F)
      seurat.obj <- regress_on_pca(seurat.obj, 'nCount_ATAC')
      seurat.obj <- FindNeighbors(seurat.atac, dims = 1:nREDUCTION, reduction = 'pca')
      seurat.obj <- FindClusters(seurat.atac, resl = 0.6)
      
}
seurat.obj <- RunUMAP(seurat.obj, reduction = "pca", dims = 1:nREDUCTION)


if(integrate_by == 'harmony'){
    library(harmony)
    seurat.obj <- RunHarmony(seurat.obj, c("sample"), assay.use = 'ATAC')

    seurat.obj <- seurat.obj %>% 
        RunUMAP(reduction = "harmony", dims = 1:nREDUCTION) 
}


## clustering by louvain algorithm
  ## seurat implemented louvain algorithm
  redm = ifelse(integrate_by == 'harmony', 'harmony', 'pca')
  seurat.obj = FindNeighbors(seurat.obj, reduction = redm, dims = 1:nREDUCTION, k.param = 20)
  if (toupper(k) == 'NULL' || k == '0'){
    resl = 0.2
  }else{
    k = as.integer(k)
    resl = queryResolution4Seurat(seurat.obj, reduction = 'pca', npc = nREDUCTION, k = k,
                                min_resl = 0.01)
  }
  seurat.obj = FindClusters(seurat.obj, resolution = resl)
  seurat.obj$active_clusters = seurat.obj$seurat_clusters


saveRDS(seurat.obj, file = paste0(output_dir, '/seurat_obj_', integrate_by, '.rds'))

#output bacrcode cluster information
bc_cls = data.table('Barcode' = rownames(seurat.obj@meta.data), 'Cluster' = seurat.obj@meta.data$active_clusters)
setkey(bc_cls, Cluster)

write.table(bc_cls, file = paste0(output_dir, '/cell_cluster_table_', integrate_by, '.txt'), sep = '\t',
            quote = F, row.names = F)

cg <- DimPlot(seurat.obj, reduction = 'umap', group.by = 'active_clusters', label = T) + 
  theme(legend.text = element_text(size = 17)) 
if(length(unique(seurat.obj$active_clusters)) < 10) cg = cg + scale_color_brewer(palette = "Paired")

cg1 <- DimPlot(seurat.obj, reduction = 'umap', group.by = 'sample') + 
  theme(legend.text = element_text(size = 17)) + scale_color_brewer(palette = "Set1")

pfname = paste0(output_dir, '/umap_clusters_', integrate_by, '.eps')
   ggsave(CombinePlots(plots = list(cg1, cg)), file = pfname, device = 'eps', width = 14, height = 6)

