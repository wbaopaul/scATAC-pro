library(bedr)
library(data.table)
options(scipen = 999)

args = commandArgs(T)

inputs = args[1]
merged_file = args[2]

tmp = unlist(strsplit(inputs, ','))
gap = as.integer(tmp[length(tmp)])
if(is.na(gap)) {
    ## gap not specified
    gap = 200
    files = tmp
}else{
    files = tmp[-length(tmp)]
}

rm(tmp)
peaks = NULL

for(file0 in files){
    peaks = rbind(peaks, fread(file0))
}

peaks = peaks[!grepl(V1, pattern = 'random', ignore.case = T)]
peaks = peaks[!grepl(V1, pattern = 'Un', ignore.case = T)]

regions = paste0(peaks$V1, ':', peaks$V2, '-', peaks$V3)
a.sort   <- bedr.sort.region(regions)
a.merged <- bedr.merge.region(a.sort, distance = gap)

dd = data.table('peak' = a.merged)
dd[, 'chr' := unlist(strsplit(peak, ':'))[1], by = peak]
dd[, 'range' := unlist(strsplit(peak, ':'))[2], by = peak]
dd[, 'start' := unlist(strsplit(range, '-'))[1], by = peak]
dd[, 'end' := unlist(strsplit(range, '-'))[2], by = peak]

dd[, c('peak', 'range') := NULL]

write.table(dd, file = merged_file, quote = F, 
    row.names = F, col.names = F, sep = '\t')


