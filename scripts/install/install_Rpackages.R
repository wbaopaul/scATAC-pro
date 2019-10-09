## install R packages for scATAC-pro

if(!require(BiocManager)){
    install.packages('BiocManager')
}

pks = c('devtools', 'flexdashboard', 'png', 'data.table', 'Matirx', 'Rcpp', 'ggplot2', 'flexmix',
  'optparse', 'magrittr', 'readr', 'Seurat', 'bedr', 'gridExtra', 'ggrepel', 'kableExtra', 'viridis')

for(pk in pks){
    if(!require(pk, character.only = T)) {
        message(paste('Install', pk, '...'))
        install.packages(pk, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    }
}

bioc.pks = c('RColorBrewer','pheatmap','motifmatchr', 'chromVAR', 'SummarizedExperiment', 'BiocParallel', 'DESeq2', 'edgeR', 'matrixStats')

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
