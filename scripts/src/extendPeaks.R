library(data.table)

args = commandArgs(T)
peakFile = args[1]

dd = fread(peakFile)

names(dd)[2:3] = c('start', 'end')
dd[, 'ss' := (end - start)]
dd[, 'midp' := floor(end/2 + start/2)]
dd[, 'start' := ifelse(ss < 500, midp - 250, start)]
dd[, 'end' := ifelse(ss < 500, midp + 250, end)]
dd[, c('ss', 'midp') := NULL]

write.table(dd, file = peakFile, row.names = F, col.names = F, quote = F,
             sep = '\t')
