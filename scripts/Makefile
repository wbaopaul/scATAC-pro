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
	@$(SCRIPTS)/dex_fastq.sh $(INPUT_FILE) $(OUTPUT_DIR)
	@echo "Demuxtiplex Done"




######################################
## bwa Alignment
##
######################################

mapping: /bin/bash
	@echo "--------------------------------------------" 
	@date 
	@echo "Aligning by bwa ..." 
	@$(SCRIPTS)/mapping_bwa.sh $(INPUT_FILE) $(OUTPUT_DIR)
	@echo "Alignment Done"


#######################################
## filter bam
#######################################
filter_bam: /bin/bash
	@echo "--------------------------------------------" 
	@date 
	@echo "Filtering bwa ..." 
	@$(SCRIPTS)/filter_bam.sh $(INPUT_FILE) $(OUTPUT_DIR)
	@echo "Filtering Bam Done"




#######################################
## call peak
#######################################
call_peak: 
	@echo "--------------------------------------------" 
	@date 
	@echo "Call peaks ..." 
	@$(SCRIPTS)/call_peak.sh $(INPUT_FILE) $(OUTPUT_DIR)
	@echo "Call peaks Done"


