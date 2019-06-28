library(zinba)

generateAlignability(
  
  mapdir='/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/Mappability_Map/map50_hg38/', #mappability directory from unpacked mappability files
  
  outdir='/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/output/peaks/ZINBA/', #directory for processed files, used later in analysis
  
  athresh=1, #number of hits per read allowed during mapping process
  
  extension=200, #average fragment library length
  
  twoBitFile='/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/hg38.2bit' #path to downloaded genome build file in .2bit format
  
)


basealigncount(
  
  inputfile='/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/output/fragments/tmp.bed',
  
  outputfile='/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/output/peaks/ZINBA/ss.bwa.basecount',
  
  extension=200,
  
  filetype='bed',
  
  twoBitFile='/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/hg38.2bit'
  
)

zinba(
  
  align='/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/output/peaks/ZINBA/',
  
  numProc=4,
  
  seq='/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/output/fragments/tmp.bed',
  
  basecountfile='/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/output/peaks/ZINBA/ss.bwa.basecount',
  
  filetype="bed",
  
  outfile="/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/output/peaks/ZINBA/",
  
  twoBit='/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/hg38.2bit',
  
  extension=200,
  
  printFullOut=1,
  
  refinepeaks=1,
  
  broad=F,
  
  input="none"
  
)



