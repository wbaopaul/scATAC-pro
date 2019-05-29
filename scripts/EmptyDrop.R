library(DropletUtils)
library(data.table)
library(Matrix)

args = commandArgs(T)

input_mtx_file = args[1]
output_pre = args[2]
fdr = as.numeric(args[3])

mat = fread(input_mtx_file)
features = mat[, 1]
mat0 = Matrix(as.matrix(mat[, -1]), sparse = T)
rm(mat)
  set.seed(2019)
  cell.out <- emptyDrops(mat0)
  
  filter.out <- cell.out[complete.cases(cell.out), ]
  
  
  filter.out = filter.out[filter.out$FDR <= fdr, ]
  
  select.cells = rownames(filter.out)

  filter.knee = filter.out[filter.out$FDR == 0, ]
  select.cells.knee = rownames(filter.knee)

out_mat = cbind(features, mat0[, colnames(mat0) %in% select.cells]) 
write.table(out_mat, file = paste0(output_pre, '_fdr', fdr), sep = '\t', row.names = F, quote = F)  

out_mat0 = cbind(features, mat0[, colnames(mat0) %in% select.cells.knee]) 
write.table(out_mat0, file = paste0(output_pre, '_kneepoint'), sep = '\t', row.names = F, quote = F)  



