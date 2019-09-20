library(data.table)

args = commandArgs(T)
fileName = args[1] ## cell cluster table
cut_N = as.integer(args[2])

dirName = dirname(fileName)
name0 = basename(fileName)
dd = fread(fileName)

dd[, 'N' := .N, by = Cluster]

dd = dd[N >= cut_N]
dd[, 'N' := NULL]

write.table(dd, file = paste0(dirName, '/filtered_', name0), sep = '\t',
   row.names = F, quote = F)
