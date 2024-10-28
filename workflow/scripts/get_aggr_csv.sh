#!/bin/bash

find $1 -iname molecule_info.h5 | \
awk 'BEGIN{print "sample_id,molecule_h5,original_sample_id"}
{
    n=split($0,s,"/"); 
    sample_id=s[n-2]; 
    gsub("_count","",sample_id); 
    print sample_id","$0","sample_id;
}'