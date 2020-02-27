library(data.table)
library(ggrepel)
library(gridExtra)
library(pheatmap)
library(viridis)

args = commandArgs(T)
group1_fp = args[1]
group2_fp = args[2]
pvalue_fp = as.numeric(args[3])
fp_res_dir = args[4]
output_dir = args[5]


process_fp_stat <- function(footprint_stats, pvalue_fp, cl1, cl2){
  footprint_stats[, 'motif1' := unlist(strsplit(Motif, '.', fixed = T))[3], by = Motif]
  footprint_stats[, 'motif2' := unlist(strsplit(Motif, '.', fixed = T))[4], by = Motif]
  footprint_stats[, 'motif2' := ifelse(is.na(motif2), "", motif2)]
  footprint_stats[, 'motif' := paste0(motif1, motif2)]
  footprint_stats[, c('motif1', 'motif2') := NULL]
  
  footprint_stats[, 'isSig' := ifelse(P_values <= pvalue_fp, 'differentiated', 'no difference')]
  footprint_stats[, 'isSig' := ifelse(P_values <= pvalue_fp & TF_Activity > 0, paste0('cluster_', cl2, '_higher'), isSig)]
  footprint_stats[, 'isSig' := ifelse(P_values <= pvalue_fp & TF_Activity < 0, paste0('cluster_', cl1, '_higher'), isSig)]
  
  footprint_stats$motif_show = ""
  footprint_stats$motif = toupper(footprint_stats$motif)
  footprint_stats[P_values <= pvalue_fp]$motif_show = footprint_stats[P_values <= pvalue_fp]$motif
  return(footprint_stats)
}


if(group1_fp != 'one'){
  
  footprint_stats.file = paste0(fp_res_dir, '/', group1_fp, '_vs_', 
                                group2_fp, '/differential_statistics.txt')
  if(file.exists(footprint_stats.file)){
   
    #  library(grid)
    footprint_stats = fread(footprint_stats.file)
    footprint_stats <- process_fp_stat(footprint_stats, pvalue_fp, group1_fp, group2_fp)
    
    
    p <- ggplot(data = footprint_stats, aes(x = motif, y = TF_Activity, 
                                            colour = factor(isSig), label = motif_show)) +
      geom_point() + xlab("") + 
      ylab('TF Activity Difference') + 
      theme(legend.text = element_text(size=18, family = "Helvetica"), 
            axis.title = element_text(size = 18, family = "Helvetica"), 
            axis.text.x = element_blank(), legend.title = element_blank(),
            legend.direction = 'horizontal', 
            plot.title = element_text(size = 15, family = "Helvetica",
                                      face = 'bold'), panel.background = element_rect(fill = "white"))  + 
      geom_text_repel(force = 10) + 
      theme(plot.margin=unit(c(0.2, 0, 0.2, 0), "cm"), legend.position = 'bottom') +
      scale_color_manual(values = c('#F8766D', '#619CFF', 'gray40')) 
    
    
    
   
    pfname = paste0(output_dir, '/footprint_tf_activity_',
                    group1_fp, '_vs_', group2_fp, '.eps')
      
    ggsave(p, filename = pfname, device = 'eps', height = 6,
             width = 6)
      
    
    footprint_out = subset(footprint_stats, 
                           P_values <= pvalue_fp,
                           select = c('motif_show', 'TF_Activity', 'Z_score', 'P_values', 'isSig'))

    names(footprint_out)[c(1, 5)] = c('motif', 'higher_in_cluster')
  }
  
}else{
  if(group2_fp != 'rest') stop('group2_fp must be rest if group1_fp is one')
  files = dir(fp_res_dir)
  files = files[grepl(files, pattern = '_vs_rest')]
  footprint_out = NULL
  for(dir0 in files){
    file0 = paste0(fp_res_dir, '/', dir0, '/differential_statistics.txt')
    cl0 = gsub('_vs_rest', '', dir0)
    if(!file.exists(file0)) next
    footprint_stats = fread(file0)
    footprint_stats <- process_fp_stat(footprint_stats, pvalue_fp, cl0, 'rest')
    footprint_out0 = subset(footprint_stats, 
                           P_values <= pvalue_fp,
                           select = c('motif_show', 'TF_Activity', 'Z_score', 'P_values', 'isSig'))
    
    names(footprint_out0)[c(1, 5)] = c('motif', 'higher_in_cluster')
    footprint_out = rbind(footprint_out, footprint_out0)
  }
  footprint_out = footprint_out[!grep(higher_in_cluster, pattern = 'rest_higher')]
}

## for large groups, randomly show 10 TFs
if(length(unique(footprint_out$motif)) > 100){
  footprint_out[, 'N' := .N, by = higher_in_cluster]
  cls = unique(footprint_out[N > 10]$higher_in_cluster)
  if(length(cls) >= 1){
    res0 = NULL
    for(cl0 in cls){
      tmp = footprint_out[higher_in_cluster == cl0]
      tmp = tmp[order(P_values)][1:10, ]
      res0 = rbind(res0, tmp)
    }
    footprint_out = rbind(footprint_out[N < 10], res0)
  }
}

mm = reshape2::acast(motif ~ higher_in_cluster, data = footprint_out, 
                     value.var = "P_values")
mm = -log10(mm)
mm[is.na(mm)] = 0
cn = colnames(mm)
cn.new = sapply(cn, function(x) gsub('_higher', '', x))
colnames(mm) = cn.new
mm[mm > 3] = 3
p1 <- pheatmap::pheatmap(mm, cluster_cols = F, fontsize = 13, fontsize_row = 9,
                   color = viridis::viridis(100))

pfname1 =paste0(output_dir, '/heatmap_differential_TF_footprint_', 
               group1_fp, '_vs_', group2_fp, '.eps')

ggsave(p1, filename = pfname1, device = 'eps', height = 10,
       width = 6)

write.table(footprint_out, 
            file = paste0(output_dir, '/differential_TF_footprint_', 
                          group1_fp, '_vs_', group2_fp, '.txt'),
            quote = F, row.names = F, col.names = T, sep = '\t')
  
