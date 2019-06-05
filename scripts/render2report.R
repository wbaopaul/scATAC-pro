
args = commandArgs(T)

output_file = args[1]
bc_stat_file = args[2]
filtered_bc_file = args[3]
curr_dir = args[4]

rmarkdown::render(paste0(curr_dir, "/barcode_qc_report.Rmd"), output_file=output_file,
                     params = list(bc_stat_file=bc_stat_file, selected_bcs=filtered_bc_file))
