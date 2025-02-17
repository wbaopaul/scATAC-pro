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
tss_path = args[5]

markers = fread(de_file)

de_basename = basename(de_file)

## do GO analysis ####
if(!require('clusterProfiler')){
    BiocManager::install('clusterProfiler')
}
library(clusterProfiler)
## genes with tss within DA for each cluster
tss_ann <- fread(tss_path, header = F)
names(tss_ann)[c(1:4)] <- c('chr', 'start', 'end', 'gene_name')


#markers = markers[grepl(peak, pattern = 'Tss')]
markers[, 'genes' := paste(unlist(strsplit(peak, ','))[-1], collapse = ','), by = peak]

## do go cluster by cluster
## genes associated with all peaks are set as background
## set background genes
tmp <- readRDS(paste0(output_dir, '/seurat_obj.rds'))
all.peaks = assignGene2Peak_coords(rownames(tmp), tss_ann)
bg_genes <- lapply(all.peaks, function(x) unlist(strsplit(x, ','))[-1])
rm(tmp)
bg_genes = do.call('c', bg_genes)
bg_genes = unique(sapply(bg_genes, function(x) gsub('-Tss', '', x)))

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
go_out_file = gsub('.tsv', '', go_out_file, fixed = T)
organism = ifelse(grepl(GENOME_NAME, pattern = 'mm'), 'mmu', 'hsa')
for(cl0 in cls){
    #bg_genes = unique(do.call('c', genesInDA))
    markers0 = markers[cluster == cl0]$genes
    genes0 = lapply(markers0, function(x) unlist(strsplit(x, ',')))
    if(length(genesInDA[[paste0('cluster', cl0)]]) <= 10) {
      tmp = data.table(matrix(1:8, 1, 8))
      names(tmp) = c('ID', 'Description', 'GeneRation', 'BgRatio', 'pvalue',
                     'p.adjust', 'qvalue', 'geneID')
      tmp = tmp[pvalue < 1]
      cl0 = gsub(':', '_', cl0, fixed = T)
      goByCl[[paste0('cluster', cl0)]] = tmp
      
    }else{
      tmp = do_GO(genesInDA[[paste0('cluster', cl0)]],
                  bg_genes = bg_genes,
                  type = GO_TYPE, qCutoff = 0.1, organism = organism)
      
      
        cl0 = gsub(':', '_', cl0, fixed = T)
        interm_goByCl[[paste0('cluster', cl0)]] = tmp
        goByCl[[paste0('cluster', cl0)]] = tmp@result[tmp@result$qvalue <= 0.1, ] 
      }
}
if(!require('writexl')){
    install.packages('writexl', dependencies = TRUE, repos = "http://cran.us.r-project.org")
}
library(writexl)
write_xlsx(goByCl, path = go_out_file)




