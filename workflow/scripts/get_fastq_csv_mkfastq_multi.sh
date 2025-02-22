#!/bin/bash

fastq_directory=$1
library_type=$2

find $fastq_directory -iname "*.gz" |
    grep -Ev "Undetermined|\_I1_|\_I2_" |
    grep -v fastq_path |
    sed 's:_R1_::g' |
    sed 's:_R2_::g' |
    sort |
    uniq |
    awk -v lt="$library_type" 'BEGIN{ print "fastq,sample,library_type";}
{
    n1=split($0,parts1,"_S[0-9+]_"); 
    n2=split(parts1[1],parts2,"/"); 
    fastq=$0;
    n3=split(fastq,parts3,"/");
    gsub(parts3[n3],"",fastq); 
    sample=parts2[n2]; 
    print fastq","sample","lt;
}'
