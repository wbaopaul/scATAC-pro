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


filter.out = filter.out[filter.out$FDR <= fdr, ]

select.cells = rownames(filter.out)

filter.knee = filter.out[filter.out$FDR == 0, ]
select.cells.knee = rownames(filter.knee)

out_mat = mat[, colnames(mat) %in% select.cells]
barcodes = colnames(out_mat)


dir1 = paste0(output_dir, '/fdr', fdr)
system(paste('mkdir -p', dir1))
writeMM(out_mat, file = paste0(dir1, '/matrix.mtx'))  
write.table(barcodes, file = paste0(dir1, '/barcodes.txt'), sep = '\t', 
            row.names = F, quote = F, col.names = F)
write.table(features, file = paste0(dir1, '/features.txt'), sep = '\t',
            row.names = F, quote = F, col.names = F)


out_mat0 = mat[, colnames(mat) %in% select.cells.knee]
barcodes0 = colnames(out_mat0)
writeMM(out_mat0, file = paste0(output_pre, 'filtered.kneepoint', '.mtx'))  
write.table(barcodes0, file = paste0(output_pre, 'filtered.kneepoint.barcodes.txt'), sep = '\t', row.names = F, quote = F, col.names = F)

dir2 = paste0(output_dir, '/kneepoint')
system(paste('mkdir -p', dir2))
writeMM(out_mat0, file = paste0(dir2, '/matrix.mtx'))  
write.table(barcodes0, file = paste0(dir2, '/barcodes.txt'), sep = '\t', 
            row.names = F, quote = F, col.names = F)
write.table(features, file = paste0(dir2, '/features.txt'), sep = '\t',
            row.names = F, quote = F, col.names = F)

