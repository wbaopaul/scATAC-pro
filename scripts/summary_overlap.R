library(data.table)

args = commandArgs(T)
qc_dir = args[1]
file_prefix = args[2]


## read overlap stat files
over.peaks = fread(paste0(qc_dir, '/', file_prefix, '.overlapWith.peaks'))
over.promoters = fread(paste0(qc_dir, '/', file_prefix, '.overlapWith.promoters'))
over.enhs = fread(paste0(qc_dir, '/', file_prefix, '.overlapWith.enhancers'))
freq.mito = fread(paste0(qc_dir, '/', file_prefix, '.freq.mito'))




