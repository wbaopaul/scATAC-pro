# correct barcodes using the fragment file

library(data.table)
library(magrittr)
library(Rcpp)
library(stringdist)
library(compiler)

## get reverse compliment of a sequence
rev.comp <- function(x, rev=TRUE){
  x <- toupper(x)
  y <- rep("N", nchar(x))
  xx <- unlist(strsplit(x, NULL))
  for (bbb in 1:nchar(x))
  {
    if(xx[bbb]=="A") y[bbb]<-"T"    
    if(xx[bbb]=="C") y[bbb]<-"G"    
    if(xx[bbb]=="G") y[bbb]<-"C"    
    if(xx[bbb]=="T") y[bbb]<-"A"
  }
  if(!rev) 
  {
    for(ccc in (1:nchar(x)))
    {
      if(ccc==1) yy <- y[ccc] else yy <- paste(yy, y[ccc], sep="")
    }
  }
  if(rev)
  {
    zz<-rep(NA, nchar(x))
    for(ccc in (1:nchar(x)))
    {
      zz[ccc] <- y[nchar(x)+1-ccc]
      if(ccc == 1) yy <- zz[ccc] else yy <- paste(yy, zz[ccc],sep="")
    }
  }
  return(yy)  
}
rev.comp = cmpfun(rev.comp)


args = commandArgs(T)

frag.file = args[1]
whitelist.file = args[2]
mis.match.cutoff = as.integer(args[3]) # 1 as default
corrected.frag.file = args[4]
doReverseComp = as.logical(args[5])

frags = fread(frag.file, header = F)
names(frags) = c('chr', 'start', 'end', 'bc', 'ndup')
white.list = fread(whitelist.file, header = F)
names(white.list) = 'wbc'

if(doReverseComp) white.list$wbc = sapply(white.list$wbc, rev.comp)

setkey(frags, bc)
setkey(white.list, wbc)

frags1 = frags[bc %in% white.list$wbc]
frags0 = frags[!bc %in% white.list$wbc]
rm(frags)

if(mis.match.cutoff == 0){
  write.table(frags1, file = corrected.frag.file, sep = '\t',
              row.names = F, quote = F)
}else{
  
  # get the unique frequence of each barcode 
  
  freq.white = subset(frags1, select = 'bc')
  freq.white[, 'N' := .N, by = bc]
  freq.white %<>% .[!duplicated(.)]
  freq.white %<>% .[order(-N), ]
  # for each non white barcode, get best matched bc in white list
  freq.nonwhite = subset(frags0,  select = 'bc')
  freq.nonwhite[, 'N' := .N, by = bc]
  freq.nonwhite %<>% .[!duplicated(.)]
  
  nbc = freq.nonwhite$bc
  
  #nbc2.wbc = barcode_match(nbc, freq.white$bc, 1)
  nbc2.wbc = sapply(nbc, function(x) stringdist::amatch(x, freq.white$bc, 
                                          maxDist=1, method = 'hamming'))
  
  match.table = data.table('nbc' = names(nbc2.wbc), 'wbc' = nbc2.wbc) 
  match.table %<>% .[!is.na(wbc)]
  match.table$wbc = freq.white$bc[match.table$wbc]
  
  frags0 %<>% .[bc  %in% match.table$nbc]
  setkey(match.table,  nbc)
  frags0[, 'wbc' := match.table[J(frags0$bc)]$wbc]
  frags0[, 'bc' := NULL]
  setcolorder(frags0, c('chr', 'start', 'end', 'bc', 'ndup'))
  frags = rbind(frags1, frags0)
  rm(frags1, frags0)
  frags[, 'tmark' := paste(chr, start, end, bc, sep = '-')]
  frags[, 'ndup1' := sum(ndup), by = tmark]
  frags[, c('ndup', 'tmark' := NULL)]
  write.table(frags, file = corrected.frag.file, sep = '\t',
              row.names = F, quote = F)
}




