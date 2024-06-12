
args = commandArgs(T)

output_file = args[1]
result_dir = args[2]
configure_file = args[3]
scatacpro_version = args[4]
sample_name = args[5]

argv <- commandArgs(trailingOnly = FALSE)
curr_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
out_dir = dirname(output_file)
library(shiny)
library(rmarkdown)
library(flexdashboard)

run(paste0(curr_dir, "/scATAC-pro_report_dynamic.Rmd"),
    shiny_args = list(launch.browser = TRUE),
    render_args = list(output_file=output_file,
                      intermediates_dir = out_dir,
                      params = list(set_title = scatacpro_version,
                                set_sample = sample_name,
                                output_dir = result_dir,
                                configure_user = configure_file)
                       )
)
