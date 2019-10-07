## scATAC-pro


## DO NOT EDIT THE REST OF THIS FILE!!

MK_PATH = $(shell dirname $(abspath $(lastword $(MAKEFILE_LIST))))

VNUM = $(shell $(MK_PATH)/scATAC-pro --version | cut -d " " -f 3)

softPATH = $(shell dirname $(shell which scATAC-pro))

INST_SCRIPTS=$(MK_PATH)/scripts
INST_SOURCES=$(INST_SCRIPTS)/src
CONFIGURE_OUT=$(wildcard ./configure_system.txt)
#CONFIG_INSTALL=$(wildcard ./configure_install.txt)
prefix ?= /usr/local/bin

install : config_check cp

######################################
## Config file
##
######################################
config_check:
ifneq ("$(CONFIGURE_OUT)","")
include $(CONFIGURE_OUT)
#include $(CONFIG_INSTALL)
else
	$(error configure_system.txt file not found. Please run 'make configure' first)
endif

######################################
## Configure
##
######################################
configure:
ifndef prefix
	$(error prefix not defined !)
else
	make -f ./scripts/install/Makefile PREFIX=$(realpath $(prefix))
endif

######################################
## Compile
##
######################################



test: config_check
	@echo ${PYTHON_PATH}

######################################
## Create installed version
##
######################################

cp:	
	@echo "Installing ..."
	$(eval INSTALL_PATH := $(realpath $(INSTALL_PREFIX))/scATAC-pro_$(VNUM))
	@if [ -d $(INSTALL_PATH) ]; then \
		rm -rf $(INSTALL_PATH); \
	fi

ifneq ($(realpath $(MK_PATH)), $(realpath $(INSTALL_PATH)))
	@cp -Ri $(MK_PATH) $(INSTALL_PATH)
endif
	@echo export PATH=$(INSTALL_PATH):"$$"PATH >> ~/.bashrc
	@echo "scATAC-pro installed in $(shell realpath $(INSTALL_PATH)) !"
	@bash


