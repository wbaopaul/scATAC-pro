
args = commandArgs(T)

output_file = args[1]
result_dir = args[2]
configure_file = args[3]

argv <- commandArgs(trailingOnly = FALSE)
curr_dir <- dirname(substring(argv[grep("--file=", argv)], 8))

rmarkdown::render(paste0(curr_dir, "/scATAC-pro_report.Rmd"), output_file=output_file,
                     params = list(output_dir = result_dir, configure_user = configure_file)) 
