
## generate gene cis activity score using unique fragment file (not from the peak-cell count matrix) 
## the total reads cutoff
library(data.table)
library(Rcpp)
library(Matrix)



#sourceCpp(paste0('getOverlaps.cpp'))

sourceCpp(code='
  #include <Rcpp.h>
  using namespace Rcpp;
  // for each read, return which region it belongs to; 0 indicates not overlaped with
  // any region
  // should sort both reads and regions
  // [[Rcpp::export]]
  IntegerVector getOverlaps_read2region(DataFrame reads, DataFrame regions) {
      NumericVector start1 = reads["start"];
      NumericVector end1 = reads["end"];
      NumericVector start2 = regions["start"];
      NumericVector end2 = regions["end"];
          
      int n1 = start1.size(), n2 = start2.size();
      NumericVector midP1(n1), len1(n1), len2(n2), midP2(n2);
      IntegerVector over1(n1);
          
      len1 = (end1 - start1)/2;
      midP1 = (end1 + start1)/2;
          
      len2 = (end2 - start2)/2;
      midP2 = (end2 + start2)/2;
      int k = 0;
      for(int i=0; i<n1; i++){
          over1[i] = 0;
          for(int j=k; j<n2; j++){
             if((fabs(midP1[i] - midP2[j]) <= (len1[i]+len2[j]))){
                over1[i] = j + 1;
                k = j;
                break;
              }
          }
      }
          
         return(over1);
    }'
)



args = commandArgs(T)
frags.file = args[1]
gtf.file = args[2]
out_dir = args[3]
cell_bc_file = args[4] 
include_body = as.logical(args[5])

frags = fread(frags.file, select=1:5, header = F)
names(frags) = c('chr', 'start', 'end', 'bc', 'ndup')
setkey(frags, chr, start)
frags[, 'ndup' := NULL]

barcodes = fread(cell_bc_file, header = F)
frags = frags[bc %in% barcodes$V1, ]

## read gtf file
gene_ann = fread(gtf.file, sep = '\t')
gene_ann = gene_ann[V3 == 'gene']
gene_ann[, 'gene_name' := unlist(strsplit(V9, ';'))[3], by = V9]
gene_ann[, 'gene_name' := gsub("\"", "", gene_name), by = gene_name]
gene_ann[, 'gene_name' := unlist(strsplit(gene_name, ' '))[3], by = gene_name]
names(gene_ann)[1] = 'chr'
gene_ann = subset(gene_ann, select = c(chr, V4, V5, V7, gene_name))
chrs = 1:22
chrs = c(chrs, 'X', 'Y')
gene_ann = gene_ann[chr %in% chrs]
gene_ann = gene_ann[!duplicated(gene_name)]

names(gene_ann)[2:4] = c('ss', 'ee', 'strand')
if(!include_body){
  gene_ann[, 'start' := ifelse(strand == '+', ss - 2000, ee - 2000)]
  gene_ann[, 'end' := ifelse(strand == '+', ss + 2000, ee + 2000)]
}else{
  gene_ann[, 'start' := ifelse(strand == '+', ss - 2000, ss)]
  gene_ann[, 'end' := ifelse(strand == '+', ee, ee + 2000)]
  
}

gene_ann[, 'chr' := paste0('chr', chr)]

regions = subset(gene_ann, select = c(chr, start, end, gene_name))

setkey(regions, chr, start)

chrs = unique(regions$chr)

fragsInRegion = NULL
for(chr0 in chrs){
  regions0 = regions[chr == chr0]
  frags0 = frags[chr == chr0]
  readsInRegion = getOverlaps_read2region(frags0, regions0)
  frags0[, 'region' := readsInRegion]
  frags0 = frags0[region != 0]
  
  frags0[, 'region' := regions0$gene_name[region]]
  fragsInRegion = rbind(fragsInRegion, frags0)
}

# get barcode region count matrix
mtx = subset(fragsInRegion, select = c('bc', 'region'))
setkey(mtx, bc, region)
mtx[, 'count' := .N, by = list(bc, region)]
mtx = mtx[!duplicated(mtx), ]
mtx = mtx[, 'n_per_region' := .N, by = region]
mtx = mtx[, 'n_per_bc' := .N, by = bc]
mtx = mtx[n_per_region > 5]
mtx[, c('n_per_region', 'n_per_bc') := NULL]

# covert to sparse matrix
features = sort(unique(mtx$region))
bcs = sort(unique(mtx$bc))
mtx$rowID = match(mtx$region, features)
mtx$colID = match(mtx$bc, bcs)

smat = sparseMatrix(i = mtx$rowID, j = mtx$colID, x = mtx$count,
                    dimnames = list(features, bcs))

if(!file.exists(out_dir)) system(paste('mkdir -p ', out_dir))
writeMM(smat, file = paste0(out_dir, '/matrix.mtx'))


write.table(rownames(smat), file = paste0(out_dir, '/features.txt'),
            row.names = F, quote = F, col.names = F)
write.table(colnames(smat), file = paste0(out_dir, '/barcodes.txt'),
            row.names = F, quote = F, col.names = F)



