library(bedr)
library(data.table)

args = commandArgs(T)

peak_cluster_dir = args[1]
chr_sizes_file = args[2]

files = dir(peak_cluster_dir)
files = files[grepl(files, pattern = "narrowPeak")]

chrs = fread(chr_sizes_file, header = F)$V1
peaks = NULL

for(file0 in files){
    peaks = rbind(peaks, fread(paste0(peak_cluster_dir, '/', file0)))
}
peaks = peaks[V1 %in% chrs]
peaks = peaks[!grepl(V1, pattern = 'random', ignore.case = T)]
peaks = peaks[!grepl(V1, pattern = 'Un', ignore.case = T)]
peaks = peaks[!grepl(V1, pattern = 'EBV', ignore.case = T)]

regions = paste0(peaks$V1, ':', peaks$V2, '-', peaks$V3)
a.sort   <- bedr.sort.region(regions)
a.merged <- bedr.merge.region(a.sort, distance = 200)

dd = data.table('peak' = a.merged)
dd[, 'chr' := unlist(strsplit(peak, ':'))[1], by = peak]
dd[, 'range' := unlist(strsplit(peak, ':'))[2], by = peak]
dd[, 'start' := unlist(strsplit(range, '-'))[1], by = peak]
dd[, 'end' := unlist(strsplit(range, '-'))[2], by = peak]

dd[, c('peak', 'range') := NULL]

write.table(dd, file = paste0(peak_cluster_dir, '/merged_peaks.bed'), quote = F, 
    row.names = F, col.names = F, sep = '\t')

