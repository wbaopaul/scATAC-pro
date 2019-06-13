library(data.table)
library(Matrix)
library(flexmix)
library(countreg)  ##install.packages("countreg", repos="http://R-Forge.R-project.org")


args = commandArgs(T)

input_mtx_file = args[1]
output_dir = args[2]
genome_size = as.numeric(args[3])
qc_per_bc_file = args[4]


## read matrix data
input_mtx_dir = dirname(input_mtx_file)
mat = readMM(input_mtx_file)

barcodes = fread(paste0(input_mtx_dir, '/barcodes.txt'), header = F)

colnames(mat) = barcodes$V1


## filter barcodes by frac_in_peak
qc_per_bc = fread(qc_per_bc_file)
peak_cov_frac = min(0.05, nrow(mat) * 1000/genome_size)
qc_sele_bc = qc_per_bc[frac_peak >= peak_cov_frac]

## substract counts due to contamination (rate 0.02)
CN = max(1, round(median(qc_sele_bc$total_frags)* 0.02))
qc_sele_bc[, 'total_frags' := total_frags -CN]
qc_sele_bc = qc_sele_bc[total_frags >= 0]

## fit two NB mixture model & using signal to noisy ratio to select cells
n_in_peak = Matrix::colSums(mat)
n_in_peak = n_in_peak[names(n_in_peak) %in% qc_sele_bc$bc]
flexmix(n_in_peak ~ 1, k = 2, model = FLXMRnegbin(theta = 1))
fm0 <- flexmix(n_in_peak ~ 1, k = 2, model = FLXMRnegbin())
prob1 = posterior(fm0)[, 1]
prob2 = posterior(fm0)[, 2]
mus = parameters(fm0)[1, ]

if(mus[1] > mus[2]){
  #odd = prob1/prob2
  odd = prob1
}else{
  #odd = prob2/prob1
  odd = prob2
}
aa = which(odd == 1)
select.cells = names(n_in_peak)[aa]
length(select.cells)

out_mat = mat[, colnames(mat) %in% select.cells]
barcodes = colnames(out_mat)
dim(out_mat)


system(paste('mkdir -p', output_dir))
writeMM(out_mat, file = paste0(output_dir, '/matrix.mtx'))  
write.table(barcodes, file = paste0(output_dir, '/barcodes.txt'), 
            sep = '\t', row.names = F, quote = F, col.names = F)
system(paste0('cp ', input_mtx_dir, '/features.txt ',  output_dir, '/'))

