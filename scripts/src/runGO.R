source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local('dsAnalysis_utilities.R')


args = commandArgs(T)
de_file = args[1]
output_dir = args[2]
configure_user = args[3]

read_conf <- function(configure_user){
  
  system(paste('grep =', configure_user, "|grep -v ^# | awk -F= '{print $1}' | awk '{$1=$1;print}' > vrs.txt "))
  
  system(paste('grep =', configure_user, "|grep -v ^# | awk -F= '{print $2}' | awk -F# '{print $1}' | awk '{$1=$1;print}' > vls.txt "))
  
  vrs = readLines('vrs.txt')
  vls = readLines('vls.txt')
  for(i in 1:length(vrs)){
    assign(vrs[i], vls[i], envir = .GlobalEnv)
  }
  system('rm vrs.txt')
  system('rm vls.txt')
}

read_conf(configure_user)


markers = fread(de_file)

## annotated genes with tss within each DA 
tss = fread(TSS)
names(tss)[1:7] = c('chr', 'start', 'end', 'gene_name', 'score', 'strand', 'gene_type')
tss = tss[gene_type %in% c('protein_coding', 'miRNA', 'lincRNA')]

markers$genes = 'No_TSS'
for(i in 1:nrow(markers)){
  tss0 = tss[chr == markers$chr[i]]
  tss0 = tss0[start >= markers$start[i] & end <= markers$end[i]]
  if(nrow(tss0) > 0) markers$genes[i] = paste(unique(tss0$gene_name), collapse = ',')
}

write.table(markers, file = paste0(de_file), sep = '\t',
            quote = F, row.names = F)


## do GO analysis ####
if(!require('clusterProfiler')){
    BiocManager::install('clusterProfiler')
}
library(clusterProfiler)
## genes with tss within DA for each cluster

markers = markers[genes != 'No_TSS']

## do go cluster by cluster, genes in other clusters as background
genesInDA = list()
cls = unique(markers$cluster)

for(cl0 in cls){
markers0 = markers[cluster == cl0]$genes
genes0 = lapply(markers0, function(x) unlist(strsplit(x, ',')))
genesInDA[[paste0('cluster', cl0)]] = unique(do.call('c', genes0))
}

goByCl = list()
go_out_file = paste0(output_dir, '/enrichedGO_by_cluster.xlsx')

organism = ifelse(grepl(GENOME_NAME, pattern = 'mm'), 'mmu', 'hsa')
for(cl0 in cls){
markers0 = markers[cluster == cl0]$genes
genes0 = lapply(markers0, function(x) unlist(strsplit(x, ',')))
goByCl[[paste0('cluster', cl0)]] = do_GO(genesInDA[[paste0('cluster', cl0)]],
                                        bg_genes = unique(do.call('c', genesInDA)),
                                        type = GO_TYPE, qCutoff = 0.05, organism = organism)


}
if(!require('writexl')){
    install.packages('writexl', dependencies = TRUE, repos = "http://cran.us.r-project.org")
}
library(writexl)
write_xlsx(goByCl, path = go_out_file)


