## MAKEFILE FOR scATAC-pro

## DO NOT EDIT THE REST OF THIS FILE!!

## special characters
comma := ,
space :=
space +=
slash := |
tick := -
undsc := _


## paralell


SCRIPTS = $(shell dirname $(abspath $(lastword $(MAKEFILE_LIST))))

## - no parallelisation

######################################
## System
##
######################################
config_check:
ifndef CONFIG_SYS
    $(error configure system not defined)
endif
ifndef CONFIG_FILE
	$(error USER CONFIG not defined)
endif

######################################
##  demultiplex
##
######################################

demplx_fastq: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Demultiplex fastq ..." 
	@$(SCRIPTS)/dex_fastq.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 


######################################
##  trimming adapter
##
######################################

trimming: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Trimming fastq ..." 
	@$(SCRIPTS)/trimming.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)


######################################
## Alignment
##
######################################

mapping: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Aligning ..." 
	@$(SCRIPTS)/mapping.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo


######################################
## mapping qc
##
######################################

mapping_qc: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Reporting mapping qc ..." 
	@$(SCRIPTS)/mapping_qc.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)


######################################
## Convert 10x positionsort bam to scATAC-pro bam
##
######################################

convert10xbam: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Aligning ..." 
	@$(SCRIPTS)/convert10xbam.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)

######################################
## process using bam file
## call_peak to end of process
######################################

process_with_bam: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "process uing bam file ..." 
	@$(SCRIPTS)/process_with_bam.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)

#######################################
## call peak
#######################################
call_peak: config_check 
	@echo "--------------------------------------------" 
	@date 
	@echo "Calling peaks ..." 
	@$(SCRIPTS)/call_peak.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 


#######################################
## recall peak
#######################################
recall_peak: config_check 
	@echo "--------------------------------------------" 
	@date 
	@echo "Recalling peaks given cell barcodes..." 
	@$(SCRIPTS)/recall_peak.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

#######################################
## get peak/bin barcode matrix
#######################################
get_mtx:  config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Getting peak/bin by barcode matrix ..." 
	@$(SCRIPTS)/get_mtx.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

#######################################
## generate signal
#######################################
aggr_signal: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Call peaks ..." 
	@$(SCRIPTS)/generate_signal.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

#######################################
## reads overlap with regions, per barcode
#######################################
qc_per_barcode: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Doing QC per barcode ..." 
	@$(SCRIPTS)/qc_per_barcode.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

#######################################
## call cell
#######################################
call_cell: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Calling cell ..." 
	@$(SCRIPTS)/call_cell.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

#######################################
## extract bam for cell barcodes
#######################################
get_bam4Cells: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Extracting bam for cell barcodes ..." 
	@$(SCRIPTS)/get_bam4Cells.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

#######################################
## process
#######################################
process: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Preprocessing ..." 
	@$(SCRIPTS)/process.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

#######################################
## process without demultiplex
#######################################
process_no_dex: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Preprocessing without demultiplexing..." 
	@$(SCRIPTS)/process_no_dex.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

#######################################
## clustering
#######################################
clustering: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Clustering ..." 
	@$(SCRIPTS)/clustering.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

#######################################
## motif analysis
#######################################
motif_analysis: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Doing motif analysis ..." 
	@$(SCRIPTS)/motif_analysis.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

#######################################
## runDA
#######################################
runDA: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Doing DA ..." 
	@$(SCRIPTS)/runDA.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

#######################################
## runGO
#######################################
runGO: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Doing GO analysis ..." 
	@$(SCRIPTS)/runGO.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

################################################
## split bam and generate track (bigwig)
################################################
split_bam: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "split bam into clusters ..." 
	@$(SCRIPTS)/split_bam2clusters.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

################################################
## footprint
################################################
footprint: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Footprint analysis ..." 
	@$(SCRIPTS)/footprint.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

#######################################
## downstream
#######################################
downstream: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Doing downstream analysis ..." 
	@$(SCRIPTS)/downstream.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

#######################################
## integrate
#######################################
integrate: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Doing integrating analysis ..." 
	@$(SCRIPTS)/integrate_peak.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

#######################################
## integrate_seu
#######################################
integrate_seu: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Doing integrating mtx by seurat ..." 
	@$(SCRIPTS)/integrate_seu.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

#######################################
## call_cnv
#######################################
call_cnv: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Calling CNV ..." 
	@$(SCRIPTS)/call_cvn.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

#######################################
## run_cicero
#######################################
runCicero: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Run cicero ..." 
	@$(SCRIPTS)/runCicero.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

#######################################
## report
#######################################
report: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Generating reports ..." 
	@$(SCRIPTS)/report.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

#######################################
## mergePeaks
#######################################
mergePeaks: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Merge peaks ..." 
	@$(SCRIPTS)/mergePeaks.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

#######################################
## reconstruct matrix
#######################################
reConstMtx: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "reconstruct peak-cell matrix ..." 
	@$(SCRIPTS)/reConstMtx.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

#######################################
## visualize
#######################################
visualize: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "interactively visualize ..." 
	@$(SCRIPTS)/visualize.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

#######################################
## all
#######################################
all: config_check
	@echo "--------------------------------------------" 
	@date 
	@echo "Run the whole pipeline ..." 
	@$(SCRIPTS)/all.sh $(INPUT_FILE) $(CONFIG_FILE) $(CONFIG_SYS)
	@echo 

