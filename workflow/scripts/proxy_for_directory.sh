#!/bin/bash

bcl_success_file=$1

for l in $(awk -F, 'NR>1{print $1}' $bcl_success_file); do touch ${l%/}.csv; done
