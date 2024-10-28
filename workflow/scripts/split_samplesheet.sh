#!/bin/bash

file1=$1
file2=$2
inputfile=$3
pattern=$4

awk -v f1="$file1" -v f2="$file2" -v p="$pattern" '
NR == 1 {
    print > f1;
    print > f2;
    next;
}
{
    line = tolower($0);
    where = match(line, p);
    if(where==0) {
        print $0 > f2;
    } else {
        print $0 > f1;
    }
}
' "$inputfile"
