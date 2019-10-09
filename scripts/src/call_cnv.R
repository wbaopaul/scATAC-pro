
set.seed(1)

args = commandArgs(T)
frag_file = args[1]
genome_name = args[2]
blacklist_file = args[3]
output_dir = args[4]

## read fragment file into grange format ####
blacklist = import.bed(blacklist_file)
frags = fread(frag_file)
names(frags) = c('chr', 'start', 'end', 'barcode', 'N')
frags.gr <- GRanges(seqnames = frags$chr, 
                    ranges = IRanges(frags$start, frags$end), 
                    barcode = frags$barcode, N = frags$N)
rm(frags)

if(grepl(genome_name, pattern = '38')){
  genomeName = 'BSgenome.Hsapiens.UCSC.hg38'
  if(!require(genomeName, character.only = T)) 
    BiocManager::install(genomeName)  
  
  genomeData = BSgenome.Hsapiens.UCSC.hg38
}


if(grepl(genome_name, pattern = '19') || grepl(genome_name, pattern = '37')){
  genomeName = 'BSgenome.Hsapiens.UCSC.hg19'
  if(!require(genomeName, character.only = T)) 
    BiocManager::install(genomeName) 
  
  genomeData = BSgenome.Hsapiens.UCSC.hg19
  
}
  

if(grepl(genome_name, pattern = 'mm9')){
  genomeName = 'BSgenome.Mmusculus.UCSC.mm9'
  if(!require(genomeName, character.only = T)) 
    BiocManager::install(genomeName) 
  
  genomeData = BSgenome.Mmusculus.UCSC.mm9
}

if(grepl(genome_name, pattern = 'mm10')){
  genomeName = 'BSgenome.Mmusculus.UCSC.mm10'
  if(!require(genomeName, character.only = T)) 
    BiocManager::install(genomeName) 
  
  genomeData = BSgenome.Mmusculus.UCSC.mm10
}

cnaObj <- call_CNV(frags.gr, genomeData, blacklist)


saveRDS(cnaObj, paste0(output_dir, "/CNV_obj.rds"))

