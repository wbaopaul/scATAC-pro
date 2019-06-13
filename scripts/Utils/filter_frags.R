## filter barcodes with very low total uniq reads, given the fragment file and 
## the total reads cutoff
library(data.table)

args = commandArgs(T)

frags.file = args[1]
total_reads_cutoff = as.integer(args[2])
out.frags.file = args[3]

frags = fread(frags.file, header = F)
names(frags) = c('chr', 'start', 'end', 'bc', 'ndup')
setkey(frags, chr, start, end)
frags[, 'total_reads_per_bc' := .N, by = bc]

frags.sele = frags[total_reads_per_bc > total_reads_cutoff]
frags.sele[, 'total_reads_per_bc' := NULL]

message(paste('Total lefted barcodes:', length(unique(frags.sele$bc))))

write.table(frags.sele, file = out.frags.file, sep = '\t',
                        row.names = F, quote = F, col.names = F)

                                         
