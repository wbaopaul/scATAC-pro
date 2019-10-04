## filter barcodes with very low total uniq reads, given the fragment file and 
## the total reads cutoff
library(data.table)
library(Rcpp)
library(Matrix)



#sourceCpp(paste0('getOverlaps.cpp'))

sourceCpp(code='
  #include <Rcpp.h>
  using namespace Rcpp;
  // for each read, return 1 if it overlap with any region, otherwise 0;
  // should sort both reads and regions
  // [[Rcpp::export]]
  IntegerVector getOverlaps_read2AnyRegion(DataFrame reads, DataFrame regions) {
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
                over1[i] = 1;
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
peaks.file = args[2]
promoters.file = args[3]
tss.file = args[4]
enhs.file = args[5]
out.frag.overlap.file = args[6]

frags = fread(frags.file, select=1:4, header = F)
names(frags) = c('chr', 'start', 'end', 'bc')
setkey(frags, chr, start)

frags[, 'total_frags' := .N, by = bc]
frags = frags[total_frags > 5]

frags = frags[!grepl(chr, pattern = 'random', ignore.case = T)]
frags = frags[!grepl(chr, pattern ='un', ignore.case = T)]

peaks = fread(peaks.file, select=1:3, header = F)
promoters = fread(promoters.file, select=1:3, header = F)
tss = fread(tss.file, select=1:3, header = F)
enhs = fread(enhs.file, select=1:3, header = F)
names(peaks) = names(promoters) = names(tss) =
  names(enhs) = c('chr', 'start', 'end')


setkey(peaks, chr, start)
setkey(promoters, chr, start)
setkey(tss, chr, start)
setkey(enhs, chr, start)

chrs = unique(frags$chr)

tss[, 'start' := start - 1000]
tss[, 'end' := end + 1000]
fragsInRegion = NULL
for(chr0 in chrs){
  peaks0 = peaks[chr == chr0]
  promoters0 = promoters[chr == chr0]
  tss0 = tss[chr == chr0]
  enhs0 = enhs[chr == chr0]
  frags0 = frags[chr == chr0]
  frags = frags[chr != chr0]
  if(nrow(peaks0) == 0){
    frags0[, 'peaks' := 0]
  }else{
    frags0[, 'peaks' := getOverlaps_read2AnyRegion(frags0, peaks0)]
  }
  
  if(nrow(promoters0) == 0){
    frags0[, 'promoters' := 0]
  }else{
    frags0[, 'promoters' := getOverlaps_read2AnyRegion(frags0, promoters0)]
  }
  
  if(nrow(tss0) == 0){
    frags0[, 'tss' := 0]
  }else{
    frags0[, 'tss' := getOverlaps_read2AnyRegion(frags0, tss0)]
  }
  
  if(nrow(enhs0) == 0){
    frags0[, 'enhs' := 0]
  }else{
    frags0[, 'enhs' := getOverlaps_read2AnyRegion(frags0, enhs0)]
  }
  
  fragsInRegion = rbind(fragsInRegion, frags0)
  message(paste(chr0, 'Done!'))
}
rm(frags)

# get barcode region count matrix
fragsInRegion[, 'isMito' := ifelse(chr == 'chrM', 1, 0)]
fragsInRegion[, c('chr', 'start', 'end') := NULL]
fragsInRegion[, 'frac_mito' := sum(isMito)/total_frags, by = bc]
fragsInRegion[, 'isMito' := NULL]
fragsInRegion[, 'frac_peak' := sum(peaks)/total_frags, by = bc]
fragsInRegion[, 'peaks' := NULL]
fragsInRegion[, 'frac_promoter' := sum(promoters)/total_frags, by = bc]
fragsInRegion[, 'promoters' := NULL]
fragsInRegion[, 'frac_tss' := sum(tss)/total_frags, by = bc]
fragsInRegion[, 'tss' := NULL]
fragsInRegion[, 'frac_enhancer' := sum(enhs)/total_frags, by = bc]
fragsInRegion[, 'enhs' := NULL]
fragsInRegion = fragsInRegion[!duplicated(fragsInRegion), ]

write.table(fragsInRegion, file = out.frag.overlap.file, sep = '\t',
            row.names = F, quote = F)



