source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local('dsAnalysis_utilities.R')


args = commandArgs(T)
de_file = args[1]
output_dir = args[2]
GENOME_NAME = args[3]
GO_TYPE = args[4]

markers = fread(de_file)

de_basename = basename(de_file)

## do GO analysis ####
if(!require('clusterProfiler')){
    BiocManager::install('clusterProfiler')
}
library(clusterProfiler)
## genes with tss within DA for each cluster

#markers = markers[grepl(peak, pattern = 'Tss')]
markers[, 'genes' := paste(unlist(strsplit(peak, ','))[-1], collapse = ','), by = peak]

## do go cluster by cluster, genes in other clusters as background
genesInDA = list()
cls = sort(unique(markers$cluster))

for(cl0 in cls){
    markers0 = markers[cluster == cl0]$genes
    genes0 = lapply(markers0, function(x) unlist(strsplit(x, ',')))
    genes0 = unique(do.call('c', genes0))
  #  genes0 = genes0[grepl(genes0, pattern = 'Tss')]
    genes0 = lapply(genes0, function(x) gsub('-Tss', '', x))
    genesInDA[[paste0('cluster', cl0)]] = unique(do.call('c', genes0))
}

goByCl = list()
interm_goByCl = list()
go_out_file = paste0(output_dir, '/enrichedGO_', de_basename, '.xlsx')
go_out_file = gsub('.txt', '', go_out_file, fixed = T)
organism = ifelse(grepl(GENOME_NAME, pattern = 'mm'), 'mmu', 'hsa')
for(cl0 in cls){
    bg_genes = unique(do.call('c', genesInDA))
    #if(length(cls) == 1) bg_genes = NULL
    markers0 = markers[cluster == cl0]$genes
    genes0 = lapply(markers0, function(x) unlist(strsplit(x, ',')))
    if(length(genesInDA[[paste0('cluster', cl0)]]) <= 10) {
      tmp = data.table(matrix(1:8, 1, 8))
      names(tmp) = c('ID', 'Description', 'GeneRation', 'BgRatio', 'pvalue',
                     'p.adjust', 'qvalue', 'geneID')
      tmp = tmp[pvalue < 1]
      goByCl[[paste0('cluster', cl0)]] = tmp
      
    }else{
      tmp = do_GO(genesInDA[[paste0('cluster', cl0)]],
                  bg_genes = bg_genes,
                  type = GO_TYPE, qCutoff = 0.1, organism = organism)
      
      
        interm_goByCl[[paste0('cluster', cl0)]] = tmp
        goByCl[[paste0('cluster', cl0)]] = tmp@result[tmp@result$qvalue <= 0.1, ] 
      }
}
if(!require('writexl')){
    install.packages('writexl', dependencies = TRUE, repos = "http://cran.us.r-project.org")
}
library(writexl)
write_xlsx(goByCl, path = go_out_file)

## save intermediate GO result in rds
#saveRDS(interm_goByCl, file = paste0(output_dir, '/enrichedGO_', de_basename, '.rds'))



