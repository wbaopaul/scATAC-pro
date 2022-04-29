library(bedr)
library(data.table)
library(GenomicRanges)
options(scipen = 999)

args = commandArgs(T)

inputs = args[1]
merged_file = args[2]

tmp = unlist(strsplit(inputs, ','))
lens = length(tmp)
gap = as.integer(tmp[lens - 1])
qscore = as.numeric(tmp[lens])

if(is.na(gap)){
    stop('Need to specific the mimum distance gap to merge peaks! \n
          Please add an integer number after the input peak file paths, separated by a comma)')
}

if(is.na(qscore)) {
    ## gap not specified
    qscore = 0.01
    files = tmp[-lens]
}else{
    files = tmp[1:(lens-2)]
}

rm(tmp)
peaks = NULL

for(file0 in files){
    peaks = rbind(peaks, fread(file0, header = F))
}
if(ncol(peaks) >= 9) peaks = peaks[V9 > -log10(qscore)]
#peaks = peaks[!grepl(V1, pattern = 'random', ignore.case = T)]
#peaks = peaks[!grepl(V1, pattern = 'Un', ignore.case = T)]

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


