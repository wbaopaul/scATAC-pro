## call cell by filtering barcodes given some qc stats cutoffs

library(optparse)
library(data.table)
library(Matrix)

## get input options
parser = OptionParser()

parser <- add_option(parser, c("-r", "--raw_mtx_file"), type="character", default='peak_barcode.mtx',
                     help="original peak barcode file [default %default]")

parser <- add_option(parser, c("-o", "--output_dir"), type="character", default='filtered',
                     help="output mtx directory [default %default]")


parser <- add_option(parser, c("-b", "--bc_stat_file"), type="character", default='NULL',
                     help="barcodes with summary stat file [default %default]")

parser <- add_option(parser, c("-u", "--min_uniq_frags"), type="integer", default=3000,
                help="minimal of total unique fragments per barcode [default %default]",
                metavar="number")

parser <- add_option(parser, c("-U", "--max_uniq_frags"), type="integer", default=50000,
                help="maximal of total unique per barcode [default %default]",
                metavar="number")


parser <- add_option(parser, c("-k", "--min_frac_peak"), type="double", default=0.05,
                help="minimal fraction of total pairs in peaks per barcode [default %default]",
                metavar="number")

parser <- add_option(parser, c("-t", "--min_frac_tss"), type="double", default=0,
                help="minimal fraction of total pairs overlapping with tss per barcode [default %default]",
                metavar="number")

parser <- add_option(parser, c("-e", "--min_frac_enhancer"), type="double", default=0,
                help="minimal fraction of total pairs in enhancer regions per barcode [default %default]",
                metavar="number")

parser <- add_option(parser, c("-p", "--min_frac_promoter"), type="double", default=0,
                help="minimal fraction of total pairs in tss per barcode [default %default]",
                metavar="number")
parser <- add_option(parser, c("-m", "--max_frac_mito"), type="double", default=0.2,
                     help="maximal fraction of total pairs in Mitocondrial per barcode [default %default]",
                     metavar="number")

opt = parse_args(parser)

## read summary qc table
mtx_file = opt$raw_mtx_file
output_dir = opt$output_dir
qc_bc_stat = fread(opt$bc_stat_file)

cut.min.frag = opt$min_uniq_frags
cut.max.frag = opt$max_uniq_frags
cut.mito = opt$max_frac_mito
cut.peak = opt$min_frac_peak
cut.tss = opt$min_frac_tss
cut.promoter = opt$min_frac_promoter
cut.enh = opt$min_frac_enhancer

qc_sele = qc_bc_stat[total_frags >= cut.min.frag & total_frags <= cut.max.frag &
                       frac_mito <= cut.mito &
                       frac_peak >= cut.peak &
                       frac_tss >= cut.tss &
                       frac_promoter >= cut.promoter &
                       frac_enhancer >= cut.enh]

mtx = readMM(mtx_file)
input_mtx_dir = dirname(mtx_file)
colnames(mtx) = fread(paste0(input_mtx_dir, '/barcodes.txt'), header = F)$V1
mtx = mtx[, colnames(mtx) %in% qc_sele$bc]

mtx.dir = dirname(mtx_file)
system(paste('mkdir -p', output_dir))
writeMM(mtx, file = paste0(output_dir, '/matrix.mtx'))
write.table(colnames(mtx), file = paste0(output_dir, '/barcodes.txt'), 
            sep = '\t', row.names = F, quote = F, col.names = F)
system(paste0('cp ', dirname(mtx_file), '/features.txt ', output_dir, '/features.txt'))


