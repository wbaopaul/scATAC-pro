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
	@echo 




######################################
## Alignment
##
######################################

mapping: 
	@echo "--------------------------------------------" 
	@date 
	@echo "Aligning ..." 
	@$(SCRIPTS)/mapping.sh $(INPUT_FILE) $(OUTPUT_DIR)
	@echo


######################################
## mapping qc
##
######################################

mapping_qc: /bin/bash
	@echo "--------------------------------------------" 
	@date 
	@echo "Reporting mapping qc ..." 
	@$(SCRIPTS)/mapping_qc.sh $(INPUT_FILE) $(OUTPUT_DIR)
#######################################
## filter bam
#######################################
filter_bam: /bin/bash
	@echo "--------------------------------------------" 
	@date 
	@echo "Filtering bwa ..." 
	@$(SCRIPTS)/filter_bam.sh $(INPUT_FILE) $(OUTPUT_DIR)
	@echo




#######################################
## call peak
#######################################
call_peak: 
	@echo "--------------------------------------------" 
	@date 
	@echo "Calling peaks ..." 
	@$(SCRIPTS)/call_peak.sh $(INPUT_FILE) $(OUTPUT_DIR)
	@echo 


#######################################
## get peak barcode matrix
#######################################
get_peak_barcode_mat:  
	@echo "--------------------------------------------" 
	@date 
	@echo "Getting peak by barcode matrix ..." 
	@$(SCRIPTS)/get_peak_barcode_mat.sh $(INPUT_FILE) $(OUTPUT_DIR)
	@echo 

#######################################
## get bin barcode matrix
#######################################
get_bin_barcode_mat:  
	@echo "--------------------------------------------" 
	@date 
	@echo "Getting bin by barcode matrix ..." 
	@$(SCRIPTS)/get_bin_barcode_mat.sh $(INPUT_FILE) $(OUTPUT_DIR)
	@echo 

#######################################
## generate signal
#######################################
generate_signal: 
	@echo "--------------------------------------------" 
	@date 
	@echo "Call peaks ..." 
	@$(SCRIPTS)/generate_signal.sh $(INPUT_FILE) $(OUTPUT_DIR)
	@echo 

#######################################
## reads overlap with regions, per barcode
#######################################
qc_per_barcode: 
	@echo "--------------------------------------------" 
	@date 
	@echo "Doing QC per barcode ..." 
	@$(SCRIPTS)/qc_per_barcode.sh $(INPUT_FILE) $(OUTPUT_DIR)
	@echo 

