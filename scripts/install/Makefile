## scATAC-pro
## DO NOT EDIT THE REST OF THIS FILE!!

MK_PATH = $(shell dirname $(abspath $(lastword $(MAKEFILE_LIST))))
VNUM = $(shell $(MK_PATH)/scATAC-pro --version | cut -d " " -f 3)


all : configure

######################################
## Config file
##
######################################
config_check:
ifndef PREFIX
		$(error PREFIX is not defined. Please run 'make configure prefix=YOUR_INSTALL_path')
else		
endif

######################################
## Dependencies
##
######################################
configure: config_check
		@$(MK_PATH)/install_dependencies.sh -o $(realpath $(PREFIX)) -q

