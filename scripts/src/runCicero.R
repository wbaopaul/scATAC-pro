source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local('dsAnalysis_utilities.R')


args = commandArgs(T)
seuratObj_file = args[1]
output_dir = args[2]
tss_file = args[3]
genome_size_file = args[4]

seurat.obj = readRDS(seuratObj_file)

tss_ann <- fread(tss_file, header = F)
#names(tss_ann)[c(1:4,7)] <- c('chr', 'start', 'end', 'gene_name', 'gene_type')
#tss_ann <- tss_ann[gene_type %in% c('miRNA', 'lincRNA', 'protein_coding'), ]
names(tss_ann)[c(1:4)] <- c('chr', 'start', 'end', 'gene_name')


seurat.obj$active_clusters = as.character(seurat.obj$active_clusters)

res = doCicero_gascore(seurat.obj, reduction = 'umap', genome_size_file, tss_ann, npc = 30)

ga_score = log1p(res$ga_score * 10000)

conns = res$conns

saveRDS(ga_score, file = paste0(output_dir, '/cicero_gene_activity.rds'))

#conns$Peak1 = assignGene2Peak_coords(conns$Peak1, tss_ann)
#conns$Peak2 = assignGene2Peak_coords(conns$Peak2, tss_ann)
write.table(conns, file = paste0(output_dir, '/cicero_interactions.txt'), row.names = F,
            sep = '\t', quote = F)


