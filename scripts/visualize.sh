#!/bin/bash


cello_dir=$1
 
# reading configure file
curr_dir=`dirname $0`
source ${curr_dir}/read_conf.sh
read_conf "$2"
read_conf "$3"

curr_dir=`dirname $0`

#${R_PATH}/Rscript -e  "source('${curr_dir}/src/launch_viscello.R'); cello_local('$cello_dir')"
${R_PATH}/Rscript -e  "library(VisCello.atac); cello('$cello_dir')"

