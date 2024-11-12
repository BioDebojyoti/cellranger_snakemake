#!/bin/bash

directory2search=$1

find $directory2search/*/outs/ -iname vdj_contig_info.pb |
    awk 'BEGIN{print "sample_id,vdj_contig_info"}
{
    n=split($0,s,"/"); 
    sample_id=s[n-2]; 
    gsub("_vdj","",sample_id); 
    print sample_id","$0;
}'
