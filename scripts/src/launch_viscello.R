library(VisCello.atac)
cello_local <- function(data_path = paste0(system.file("app", package='VisCello.atac'),"/data"), RStudio = F) {
    Sys.setenv("R_MAX_NUM_DLLS"=180)
    cat("Launching VisCello...")
    if(RStudio) {
        options(shiny.launch.browser = .rs.invokeShinyWindowViewer)
    }
    
    if(exists("eset", env = .GlobalEnv)) rm(eset, envir = .GlobalEnv)
    if(exists("clist", env = .GlobalEnv)) rm(clist, envir = .GlobalEnv)
    if(exists("r_data", env = .GlobalEnv)) rm(r_data, envir = .GlobalEnv)
    if(exists("global_config", env = .GlobalEnv)) rm(global_config, envir = .GlobalEnv)
    cat(paste0("Data folder: ",data_path))
    tryCatch({
        .GlobalEnv$global_config <- config::get(file = paste0(data_path, "/config.yml"), use_parent = F)
    }, error = function(x){
        stop("Cannot find config.yml in data folder, please check if data_path is correct.")
    })
    
    .GlobalEnv$mainTitle = paste0("VisCello - ",global_config$study_name)
    .GlobalEnv$organism = global_config$organism
    .GlobalEnv$study_info <- global_config$study_description
    .GlobalEnv$name_col = global_config$feature_name_column
    .GlobalEnv$id_col = global_config$feature_id_column
    tryCatch({
        .GlobalEnv$eset <- readRDS(paste0(data_path, "/eset.rds"))
        .GlobalEnv$clist <-readRDS(paste0(data_path, "/clist.rds"))
    }, error = function(x){
        stop("Cannot find eset or clist file in data folder, if data_path is correct.")
    })
    #shiny::runApp(system.file("app", package='VisCello.atac'), port = 6456, launch.browser = TRUE)
    shiny::runApp(system.file("app", package='VisCello.atac'), launch.browser = TRUE)
}
args = commadArgs(T)
cello_dir = args[1]
cello_local(cello_dir)
