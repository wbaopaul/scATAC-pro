library(data.table)
args = commandArgs(T)
frags.file = args[1]
message("Fast sorting by R using data.table...")
dd = fread(frags.file, header = F)
setkey(dd, V1, V2)

fwrite(dd, file = frags.file, scipen=100, col.names = F, sep = '\t')

