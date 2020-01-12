## filter barcodes with very low total uniq reads, given the fragment file and 
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
region.file = args[2]
bc.file = args[3]
out_dir = args[4]

bcs = fread(bc.file, header = F)

frags = fread(frags.file, select=1:5, header = F)
names(frags) = c('chr', 'start', 'end', 'bc', 'ndup')
setkey(frags, chr, start)
frags = frags[bc %in% bcs$V1]
frags[, 'ndup' := NULL]

rm(bcs)

regions = fread(region.file, select=1:3, header = F)
names(regions) = c('chr', 'start', 'end')
setkey(regions, chr, start)
regions = regions[!grepl(chr, pattern = 'random', ignore.case = T), ]
regions = regions[!grepl(chr, pattern = 'Un', ignore.case = T), ]

chrs = unique(regions$chr)

fragsInRegion = NULL
for(chr0 in chrs){
  regions0 = regions[chr == chr0]
  frags0 = frags[chr == chr0]
  readsInRegion = getOverlaps_read2region(frags0, regions0)
  frags0[, 'region' := readsInRegion]
  frags0 = frags0[region != 0]
  
  regions0[, 'rname' := paste(chr, start, end, sep = '-')]
  frags0[, 'region' := regions0$rname[region]]
  fragsInRegion = rbind(fragsInRegion, frags0)
}

# get barcode region count matrix
mtx = subset(fragsInRegion, select = c('bc', 'region'))
setkey(mtx, bc, region)
mtx[, 'count' := .N, by = list(bc, region)]
mtx = mtx[!duplicated(mtx), ]

# covert to sparse matrix
features = sort(unique(mtx$region))
bcs = sort(unique(mtx$bc))
mtx$rowID = match(mtx$region, features)
mtx$colID = match(mtx$bc, bcs)

smat = sparseMatrix(i = mtx$rowID, j = mtx$colID, x = mtx$count,
                    dimnames = list(features, bcs))

writeMM(smat, file = paste0(out_dir, '/matrix.mtx'))


write.table(rownames(smat), file = paste0(out_dir, '/features.txt'),
            row.names = F, quote = F, col.names = F)
write.table(colnames(smat), file = paste0(out_dir, '/barcodes.txt'),
            row.names = F, quote = F, col.names = F)



