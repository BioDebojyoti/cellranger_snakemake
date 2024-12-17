#!/bin/bash

dir=$1

find "$dir" -maxdepth 4 -type d -path "*/outs/per_sample_outs/*" |
    awk 'BEGIN{print "sample_id,sample_outs"}
{
    n=split($0,s,"/"); 
    donor_id=s[n]; 
    print donor_id","$0;
}'
