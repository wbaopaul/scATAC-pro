library(bedr)
library(data.table)
library(GenomicRanges)
options(scipen = 999)

args = commandArgs(T)

inputs = args[1]
merged_file = args[2]

tmp = unlist(strsplit(inputs, ','))
lens = length(tmp)

if(lens != 3){
    stop('It seems a older version of configure_user.txt file is using. \n
          Please specify parameter "mergePeaksWithin" and "filterPeaksQvalue" in the configure_user.txt file')
}

gap = as.integer(tmp[2])
qscore = as.numeric(tmp[3])

csv_contents = fread(tmp[1])
files = csv_contents$peak_file

rm(tmp)
peaks = NULL

for(file0 in files){
    peaks = rbind(peaks, fread(file0, header = F))
}
if(ncol(peaks) >= 9) peaks = peaks[V9 > -log10(qscore)]

names(peaks)[1:3] = c('chr', 'start', 'end')
schrs = standardChromosomes(makeGRangesFromDataFrame(peaks))
peaks = peaks[chr %in% schrs]

regions = paste0(peaks$chr, ':', peaks$start, '-', peaks$end)
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


