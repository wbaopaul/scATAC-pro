library(data.table)

args = commandArgs(T)

frags.file = args[1]
barcodes.file = args[2]
out.frags.file = args[3]

frags = fread(frags.file,  header = F)
names(frags)[4] = 'barcode'
barcodes = fread(barcodes.file, header = F)

frags.sele = frags[barcode %in% barcodes$V1]

write.table(frags.sele, file = out.frags.file, sep = '\t',
            row.names = F, quote = F, col.names = F)
