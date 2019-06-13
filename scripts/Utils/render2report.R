
args = commandArgs(T)

output_file = args[1]
mapping_qc_file = args[2]
bc_stat_file = args[3]
filtered_bc_file = args[4]
fragments_file = args[5]
downstream_dir = args[6]
tss_escore_file = args[7]

#curr_dir = args[6]

argv <- commandArgs(trailingOnly = FALSE)
curr_dir <- dirname(substring(argv[grep("--file=", argv)], 8))

rmarkdown::render(paste0(curr_dir, "/scATAC-pro_report.Rmd"), output_file=output_file,
                     params = list(bc_stat_file = bc_stat_file, 
                                   selected_bcs = filtered_bc_file,
                                   mapping_qc_file = mapping_qc_file,
                                   fragments_file = fragments_file,
                                   downstream_dir = downstream_dir,
                                   tss_escore_file = tss_escore_file))
