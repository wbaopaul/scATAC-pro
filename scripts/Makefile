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


######################################
##  demultiplex
##
######################################

demplx_fastq: /bin/bash
	@echo "--------------------------------------------" 
	@date 
	@echo "Demultiplex fastq ..." 
	@$(SCRIPTS)/dex_fastq.sh $(INPUT_FILE) $(CONFIG_FILE)
	@echo 


######################################
##  trimming adapter
##
######################################

trimming: /bin/bash
	@echo "--------------------------------------------" 
	@date 
	@echo "Trimming fastq ..." 
	@$(SCRIPTS)/trimming.sh $(INPUT_FILE) $(CONFIG_FILE)


######################################
## Alignment
##
######################################

mapping: 
	@echo "--------------------------------------------" 
	@date 
	@echo "Aligning ..." 
	@$(SCRIPTS)/mapping.sh $(INPUT_FILE) $(CONFIG_FILE)
	@echo


######################################
## mapping qc
##
######################################

mapping_qc: /bin/bash
	@echo "--------------------------------------------" 
	@date 
	@echo "Reporting mapping qc ..." 
	@$(SCRIPTS)/mapping_qc.sh $(INPUT_FILE) $(CONFIG_FILE)
#######################################
## filter bam
#######################################
filter_bam: /bin/bash
	@echo "--------------------------------------------" 
	@date 
	@echo "Filtering bwa ..." 
	@$(SCRIPTS)/filter_bam.sh $(INPUT_FILE) $(CONFIG_FILE)
	@echo




#######################################
## call peak
#######################################
call_peak: 
	@echo "--------------------------------------------" 
	@date 
	@echo "Calling peaks ..." 
	@$(SCRIPTS)/call_peak.sh $(INPUT_FILE) $(CONFIG_FILE)
	@echo 


#######################################
## get peak/bin barcode matrix
#######################################
get_mtx:  
	@echo "--------------------------------------------" 
	@date 
	@echo "Getting peak by barcode matrix ..." 
	@$(SCRIPTS)/get_mtx.sh $(INPUT_FILE) $(CONFIG_FILE)
	@echo 

#######################################
## generate signal
#######################################
generate_signal: 
	@echo "--------------------------------------------" 
	@date 
	@echo "Call peaks ..." 
	@$(SCRIPTS)/generate_signal.sh $(INPUT_FILE) $(CONFIG_FILE)
	@echo 

#######################################
## reads overlap with regions, per barcode
#######################################
qc_per_barcode: 
	@echo "--------------------------------------------" 
	@date 
	@echo "Doing QC per barcode ..." 
	@$(SCRIPTS)/qc_per_barcode.sh $(INPUT_FILE) $(CONFIG_FILE)
	@echo 

#######################################
## call cell
#######################################
call_cell: 
	@echo "--------------------------------------------" 
	@date 
	@echo "Calling cell ..." 
	@$(SCRIPTS)/call_cell.sh $(INPUT_FILE) $(CONFIG_FILE)
	@echo 

#######################################
## clustering
#######################################
clustering: 
	@echo "--------------------------------------------" 
	@date 
	@echo "Clustering ..." 
	@$(SCRIPTS)/clustering.sh $(INPUT_FILE) $(CONFIG_FILE)
	@echo 

#######################################
## motif analysis
#######################################
motif_analysis: 
	@echo "--------------------------------------------" 
	@date 
	@echo "Doing motif analysis ..." 
	@$(SCRIPTS)/motif_analysis.sh $(INPUT_FILE) $(CONFIG_FILE)
	@echo 

################################################
## split bam and generate track (bigwig)
################################################
split_bam: 
	@echo "--------------------------------------------" 
	@date 
	@echo "split bam into clusters ..." 
	@$(SCRIPTS)/split_bam2clusters.sh $(INPUT_FILE) $(CONFIG_FILE)
	@echo 

################################################
## footprint
################################################
footprint: 
	@echo "--------------------------------------------" 
	@date 
	@echo "Footprint analysis ..." 
	@$(SCRIPTS)/footprint_by_cluster.sh $(INPUT_FILE) $(CONFIG_FILE)
	@echo 

#######################################
## report
#######################################
report: 
	@echo "--------------------------------------------" 
	@date 
	@echo "Getting reports ..." 
	@$(SCRIPTS)/report.sh $(INPUT_FILE) $(CONFIG_FILE)
	@echo 

