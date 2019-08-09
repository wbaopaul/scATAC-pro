source_local <- function(fname){
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

source_local('dsAnalysis_utilities.R')


args = commandArgs(T)
output_dir = args[1]
configure_user = args[2]

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

if(file.exists(paste0(output_dir, '/seurat_obj.rds'))){
  seurat.obj = readRDS(paste0(output_dir, '/seurat_obj.rds'))
}else{
  stop('Should do clustering analysis first!')
}

confVar = 'nCount_ATAC'
if(test_use == 'wilxon' || test_use == 'DESeq2') confVar = NULL

if(group2 == 'others' || group2 == 'all') group2 = NULL
if(group1 == 'all') {
   cls = unique(seurat.obj$active_clusters)
   markers = NULL
   for(cluster0 in cls){
        mm = FindMarkers(seurat.obj, group.by = "active_clusters", ident.1 = cluster0, ident.2 = group2,
                test.use = test_use, max.cell.per.ident = 500, only.pos = T, latent.vars = confVar)
        mm$cluster = cluster0
        
        mm$peak = rownames(mm)
        markers = rbind(markers, mm)

   }
}else{
    if(is.null(group2)){
        markers = FindMarkers(seurat.obj, group.by = "active_clusters", ident.1 = group1, ident.2 = group2, 
                          test.use = test_use, max.cell.per.ident = 500, only.pos = T, latent.vars = confVar)
        markers$cluster = group1
    }else{
        markers = FindMarkers(seurat.obj, group.by = "active_clusters", ident.1 = group1, ident.2 = group2, 
                  test.use = test_use, max.cell.per.ident = 500, only.pos = F, latent.vars = confVar)
        markers$cluster = ifelse(markders$avg_logFC > 0, group1, group2)
    }
  
  markers$peak = rownames(markers)
}



## annotated genes with tss within each DA 
tss = fread(TSS)
names(tss)[1:7] = c('chr', 'start', 'end', 'gene_name', 'score', 'strand', 'gene_type')
tss = tss[gene_type %in% c('protein_coding', 'miRNA', 'lincRNA')]

markers[, 'chr' := unlist(strsplit(peak, '-'))[1], by = peak]
markers[, 'start' := as.integer(unlist(strsplit(peak, '-'))[2]), by = peak]
markers[, 'end' := as.integer(unlist(strsplit(peak, '-'))[3]), by = peak]

markers$genes = 'No_TSS'
for(i in 1:nrow(markers)){
tss0 = tss[chr == markers$chr[i]]
tss0 = tss0[start >= markers$start[i] & end <= markers$end[i]]
if(nrow(tss0) > 0) markers$genes[i] = paste(unique(tss0$gene_name), collapse = ',')
}

write.table(markers, file = paste0(output_dir, '/differential_peak_cluster_table.txt'), sep = '\t',
            quote = F, row.names = F)
