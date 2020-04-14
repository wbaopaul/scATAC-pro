## install R packages for scATAC-pro

if(!require(BiocManager)){
    install.packages('BiocManager')
}

pks = c('devtools', 'flexdashboard', 'png', 'data.table', 'Matirx', 'Rcpp', 'ggplot2', 'flexmix',
  'optparse', 'magrittr', 'readr', 'Seurat', 'bedr', 'gridExtra', 'ggrepel', 'kableExtra', 'viridis', 'writexl', 'xlsx', 'mefa4')

for(pk in pks){
    if(!require(pk, character.only = T)) {
        message(paste('Install', pk, '...'))
        install.packages(pk, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    }
}

bioc.pks = c('RColorBrewer','pheatmap','motifmatchr', 'chromVAR', 'SummarizedExperiment', 'BiocParallel', 'DESeq2', 'edgeR', 'matrixStats', 'cicero', 'farver', 'BSgenome.Hsapiens.UCSC.hg38', 'BSgenome.Mmusculus.UCSC.mm10',  'clusterProfiler')

for(pk in bioc.pks){
    if(!require(pk, character.only = T)) {  
        message(paste('Install', pk, '...'))
        BiocManager::install(pk)
    }
}

if(!require(chromVARmotifs)) {  
    message(paste('Install chromVARmotifs...'))
    devtools::install_github("GreenleafLab/chromVARmotifs", upgrade = 'never', force = T)
}

if(!require(harmony)) {  
    message(paste('Install harmony...'))
    devtools::install_github("immunogenomics/harmony", upgrade = 'never', force = T)
}

if(!require(VisCello.atac)) {
    devtools::install_github("qinzhu/VisCello", ref="VisCello-atac") # install
}
