#!/bin/bash

## given SampleSheet.csv file (since version 1.5.1) or peaks files(older version)
set -e

input_files=$1  ## sample sheet or two or more input features files for two or more samples, seperated by ,

# reading configure file
curr_dir=`dirname $0`

if [[ "${input_files}" == *.csv  ]]; then
   ${curr_dir}/integrate_csv.sh $input_files $2 $3 ## newer version 
else
   ${curr_dir}/integrate_peak.sh $input_files $2 $3 ## <= v1.5.0
fi
