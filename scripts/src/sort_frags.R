library(data.table)
args = commandArgs(T)
frags.file = args[1]

dd = fread(frags.file, header = F)
head(dd)

fwrite(dd, file = frags.file, scipen=100, col.names = F, sep = '\t')

