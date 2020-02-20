library(data.table)

args = commandArgs(T)
peak_file = args[1]
chr_sizes_file = args[2]

dd = fread(peak_file, header = F)
chrs = fread(chr_sizes_file, header = F)$V1

dd = dd[V1 %in% chrs]


write.table(dd, file = peak_file, sep = '\t', row.names = F, col.names = F, quote = F)
