## install R packages for scATAC-pro

if(!require(BiocManager)){
    install.packages('BiocManager', dependencies = TRUE, repos = "http://mran.revolutionanalytics.com/snapshot/2019-11-01")
}

pks = c('devtools', 'flexdashboard', 'png', 'data.table', 'Matirx', 'Rcpp', 'ggplot2', 'flexmix',
  'optparse', 'magrittr', 'readr', 'Seurat', 'bedr', 'gridExtra', 'ggrepel', 'kableExtra', 'viridis', 'writexl')

for(pk in pks){
    if(!require(pk, character.only = T)) {
        message(paste('Install', pk, '...'))
        install.packages(pk, dependencies = TRUE, repos = "http://mran.revolutionanalytics.com/snapshot/2019-11-01")
    }
}

bioc.pks = c('RColorBrewer','pheatmap','motifmatchr', 'chromVAR', 'SummarizedExperiment', 'BiocParallel', 'DESeq2', 'edgeR', 'matrixStats', 'cicero', 'farver', 'org.Mm.eg.db', 'org.hs.eg.db', 'BSgenome.Hsapiens.UCSC.hg38', 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Mmusculus.UCSC.mm10', 'BSgenome.Mmusculus.UCSC.mm9', 'clusterProfiler')

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

if(!require(cisTopic)) {  
    message(paste('Install cisTopic...'))
    devtools::install_github("aertslab/cisTopic", upgrade = 'never', force = T)
}
