#!/bin/bash

pat="Saving pipestance info to \"";

curr_pipestance_dir=$(awk -v p="$pat" '{if($0 ~ p) {gsub(".mri.tgz\"","",$0); n=split($0,dir,"/"); print dir[n];}}' $1);
output_directory=$2;


if [[ $output_directory != $curr_pipestance_dir ]]; then
    cp \-R $curr_pipestance_dir/* $output_directory/;        
    rm \-r $curr_pipestance_dir;
fi
