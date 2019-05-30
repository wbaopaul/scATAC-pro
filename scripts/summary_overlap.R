library(data.table)

args = commandArgs(T)
qc_dir = args[1]
file_prefix = args[2]


## read overlap stat files
over.peaks = fread(paste0(qc_dir, '/', file_prefix, '.overlapWith.peaks'))
over.promoters = fread(paste0(qc_dir, '/', file_prefix, '.overlapWith.promoters'))
over.enhs = fread(paste0(qc_dir, '/', file_prefix, '.overlapWith.enhancers'))
over.tss = fread(paste0(qc_dir, '/', file_prefix, '.overlapWith.tss'))
freq.mito = fread(paste0(qc_dir, '/', file_prefix, '.freq.mito'))

merged.dd = cbind(over.peaks, freq.mito[, 3], over.promoters[, 3], over.tss[, 3],
                  over.enhs[, 3])
names(merged.dd)[3:7] = c('read_in_peak', 'read_in_mito', 'read_in_promoter', 'read_in_tss',
                          'read_in_enhancers') 
write.table(merged.dd, file = paste0(qc_dir, '/', file_prefix, '.qc_per_bc'), quote = F,
            row.names = F, sep = '\t')
