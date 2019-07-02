## scATAC-pro


## DO NOT EDIT THE REST OF THIS FILE!!

MK_PATH = $(shell dirname $(abspath $(lastword $(MAKEFILE_LIST))))

INST_SCRIPTS=$(MK_PATH)/scripts
INST_SOURCES=$(INST_SCRIPTS)/src
CONFIGURE_OUT=$(wildcard ./configure_system.txt)
CONFIG_SYS=$(wildcard ./configure_install.txt)


install : config_check mapbuilder readstrimming iced cp

######################################
## Config file
##
######################################
config_check:
ifneq ("$(CONFIGURE_OUT)","")
include $(CONFIGURE_OUT)
else
	$(error configure_system.txt file not found. Please run 'make configure' first)
endif

######################################
## Configure
##
######################################
configure:
ifneq ("$(CONFIG_SYS)","")
	make -f ./scripts/install/Makefile CONFIG_SYS=$(CONFIG_SYS)
else
	$(error configure_install.txt file not found !)
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
ifneq ($(realpath $(MK_PATH)), $(realpath $(INSTALL_PATH)))
	cp -Ri $(MK_PATH) $(INSTALL_PATH)
endif
@echo "scATAC-pro installed in $(shell realpath $(INSTALL_PATH)) !"
