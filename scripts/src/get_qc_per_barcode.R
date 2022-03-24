## filter barcodes with very low total uniq reads, given the fragment file and 
## the total reads cutoff
library(data.table)
library(Rcpp)
library(Matrix)
library(GenomicRanges)


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
          
      len1 = (end1 - start1 + 1)/2;
      midP1 = (end1 + start1)/2;
          
      len2 = (end2 - start2 + 1)/2;
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
  }
    
// [[Rcpp::export]]
  NumericVector getOverlaps_read2AllRegion(DataFrame reads, DataFrame regions) {
      NumericVector start1 = reads["start"];
      NumericVector end1 = reads["end"];
      
      NumericVector start2 = regions["start"];
      NumericVector end2 = regions["end"];
          
      int n1 = start1.size(), n2 = start2.size();
      NumericVector midP1(n1), len1(n1), len2(n2), midP2(n2);
      NumericVector over1(n1);
          
      len1 = (end1 - start1 + 1)/2;
      midP1 = (end1 + start1)/2;
          
      len2 = (end2 - start2 + 1)/2;
      midP2 = (end2 + start2)/2;
      int j = 0;
      int k = 0;
      //for(int i=0; i<n1; i++){
      //    over1[i] = 0;
      //}
      //over1 = 0;
      if(start1[0] > end2[n2-1] || end1[n1-1] < start2[0]) {
          return(over1);
      }
      for(int i=0; i<n1; i++){
          while (k<n2 && fabs(midP1[i] - midP2[k]) > (len1[i]+len2[k])){
            k++;
          }
          // the kth element is the first element overlapped with the current fragment
          // take advantage of sorted fragments and regions, for next fragment,
          // it will start searching from kth region
          j=k;
          
          // current frag not overlapped any region then search from beginning
          if(j >= n2 -1) k = 0, j = 0;  
          while (j<n2 && fabs(midP1[i] - midP2[j]) <= (len1[i]+len2[j])){
                over1[i] = over1[i] + 1;
                j++;
          }
      }
          
         return(over1);
    }

  
 // [[Rcpp::export]]
  NumericVector getOverlaps_tss2Reads(DataFrame regions, DataFrame left_flank, DataFrame reads) {
      NumericVector start1 = regions["start"];
          
      int n = start1.size();
      NumericVector over(n), left_over(n);
      NumericVector tss_enrich_score(n);
      over = getOverlaps_read2AllRegion(regions, reads);
      left_over = getOverlaps_read2AllRegion(left_flank, reads);
      tss_enrich_score = (over+1.0) / (left_over+1.0) ;
      return(tss_enrich_score);
    }    
  
  
 '
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

## only keep reads in standard chrs
chrs = standardChromosomes(makeGRangesFromDataFrame(frags))
frags = frags[chr %in% chrs]
frags[, 'total_frags' := .N, by = bc]
frags = frags[total_frags > 5]

#frags = frags[!grepl(chr, pattern = 'random', ignore.case = T)]
#frags = frags[!grepl(chr, pattern ='un', ignore.case = T)]

peaks = fread(peaks.file, select=1:3, header = F)
tss = fread(tss.file, select=1:3, header = F)
names(peaks) = names(tss) = c('chr', 'start', 'end')

## promoters.file and enhs.file are optional
if(file.exists(promoters.file)) {
    promoters = fread(promoters.file, select=1:3, header = F)
    names(promoters) =  c('chr', 'start', 'end')
    setkey(promoters, chr, start)
}
if(file.exists(enhs.file)) {
    enhs = fread(enhs.file, select=1:3, header = F)
    names(enhs) =  c('chr', 'start', 'end')
    setkey(enhs, chr, start)
}

setkey(peaks, chr, start)
setkey(tss, chr, start)

#chrs = unique(frags$chr)

## calculate tss enrichment score
if(T){
  set.seed(2021)
  tss4escore = copy(tss)
  if(nrow(tss) > 40000) tss4escore = tss[sort(sample(1:nrow(tss), 40000))]
  setkey(tss4escore, chr, start)
  
  tss4escore = unique(tss4escore)
  tss4escore.left = copy(tss4escore)
  tss4escore[, 'start' := start - 50]
  tss4escore[, 'end' := start + 100]
  tss4escore.left[, 'start' := start - 1950]
  tss4escore.left[, 'end' := start + 100]
  
  bcs = unique(frags$bc)
  escores = rep(1, length(bcs))
  names(escores) = bcs
  frags4escore = frags[total_frags >= 1000] ## only caculate tss_escore for selected bc
  
  for(bc0 in unique(frags4escore$bc)){
    frags0 = frags4escore[bc == bc0]
    chrs0 = unique(frags0$chr)
    tss4escore0 = tss4escore[chr %in% chrs0]
    tss4escore0.left = tss4escore.left[chr %in% chrs0]
    #chrs = unique(frags0$chr)
    escores_chrs = NULL
    for(chr0 in chrs0){
      if(nrow(frags0[chr==chr0]) <= 50) next
      escores_chrs[chr0] = max(getOverlaps_tss2Reads(tss4escore0[chr==chr0], 
                                                     tss4escore0.left[chr==chr0],
                                                     frags0[chr==chr0]))
    }
    
    escores[bc0] = max(escores_chrs)
  }
  escores = escores + runif(length(escores), 0, 0.1)
  rm(frags4escore, frags0, tss4escore, tss4escore.left)
}


tss[, 'start' := start - 1000]
tss[, 'end' := end + 1000]
fragsInRegion = NULL
for(chr0 in chrs){
  peaks0 = peaks[chr == chr0]

  tss0 = tss[chr == chr0]
  frags0 = frags[chr == chr0]
  frags = frags[chr != chr0]
  if(nrow(peaks0) == 0){
    frags0[, 'peaks' := 0]
  }else{
    frags0[, 'peaks' := getOverlaps_read2AnyRegion(frags0, peaks0)]
  }
  
  if(nrow(tss0) == 0){
    frags0[, 'tss' := 0]
  }else{
    frags0[, 'tss' := getOverlaps_read2AnyRegion(frags0, tss0)]
  }

  if(file.exists(promoters.file)) {
    promoters0 = promoters[chr == chr0]
    if(nrow(promoters0) == 0){
        frags0[, 'promoters' := 0]
    }else{
        frags0[, 'promoters' := getOverlaps_read2AnyRegion(frags0, promoters0)]
    }
  }
  
  
  if(file.exists(enhs.file)){ 
    enhs0 = enhs[chr == chr0]
    if(nrow(enhs0) == 0){
        frags0[, 'enhs' := 0]
    }else{
        frags0[, 'enhs' := getOverlaps_read2AnyRegion(frags0, enhs0)]
    }
  } 
 
  
  fragsInRegion = rbind(fragsInRegion, frags0)
  message(paste(chr0, 'Done!'))
}
rm(frags)

# get barcode region count matrix
fragsInRegion[, 'isMito' := ifelse(chr %in% c('chrM'), 1, 0)]
fragsInRegion[, c('chr', 'start', 'end') := NULL]
fragsInRegion[, 'frac_mito' := sum(isMito)/total_frags, by = bc]
fragsInRegion[, 'isMito' := NULL]
fragsInRegion[, 'frac_peak' := sum(peaks)/total_frags, by = bc]
fragsInRegion[, 'peaks' := NULL]
fragsInRegion[, 'frac_tss' := sum(tss)/total_frags, by = bc]
fragsInRegion[, 'tss' := NULL]

if(file.exists(promoters.file)) {
    fragsInRegion[, 'frac_promoter' := sum(promoters)/total_frags, by = bc]
    fragsInRegion[, 'promoters' := NULL]
}else{
    fragsInRegion[, 'frac_promoter' := 0]
}

if(file.exists(enhs.file)) {
    fragsInRegion[, 'frac_enhancer' := sum(enhs)/total_frags, by = bc]
    fragsInRegion[, 'enhs' := NULL]
}else{
    fragsInRegion[, 'frac_enhancer' := 0]
}
fragsInRegion = unique(fragsInRegion)
fragsInRegion$tss_enrich_score = escores[fragsInRegion$bc]

write.table(fragsInRegion, file = out.frag.overlap.file, sep = '\t',
            row.names = F, quote = F)



