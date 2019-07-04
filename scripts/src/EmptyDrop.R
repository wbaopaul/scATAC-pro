if(!require("DropletUtils")) BiocManager::install('DropletUtils')
library(DropletUtils)
library(data.table)
library(Matrix)

args = commandArgs(T)

input_mtx_file = args[1]
output_dir = args[2]
fdr = as.numeric(args[3])

input_mtx_dir = dirname(input_mtx_file)

mat = readMM(input_mtx_file)
features = fread(paste0(input_mtx_dir, '/features.txt'), header = F)
barcodes = fread(paste0(input_mtx_dir, '/barcodes.txt'), header = F)
rownames(mat) = features$V1
colnames(mat) = barcodes$V1

set.seed(2019)
cell.out <- emptyDrops(mat)

filter.out <- cell.out[complete.cases(cell.out), ]

saveRDS(filter.out, file = paste0(output_dir, '/EmptyDrop_obj.rds'))

filter.out = filter.out[filter.out$FDR <= fdr, ]

select.cells = rownames(filter.out)

out_mat = mat[, colnames(mat) %in% select.cells]
barcodes = colnames(out_mat)


system(paste('mkdir -p', output_dir))
writeMM(out_mat, file = paste0(output_dir, '/matrix.mtx'))  
write.table(barcodes, file = paste0(output_dir, '/barcodes.txt'), sep = '\t', 
            row.names = F, quote = F, col.names = F)
write.table(features, file = paste0(output_dir, '/features.txt'), sep = '\t',
            row.names = F, quote = F, col.names = F)


