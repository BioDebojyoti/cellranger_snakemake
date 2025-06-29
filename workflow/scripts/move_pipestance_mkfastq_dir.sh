#!/bin/bash

pat="Saving pipestance info to \""

curr_pipestance_dir=$(awk -v p="$pat" '{if($0 ~ p) {gsub(".mri.tgz\"","",$0); n=split($0,dir,"/"); print dir[n];}}' $1)
output_directory=$2

if [[ $output_directory != $curr_pipestance_dir ]]; then
    cp \-r $curr_pipestance_dir $output_directory/
    rm \-r $curr_pipestance_dir
fi

fastq_files=$(find "${output_directory}"/"${curr_pipestance_dir}" -iname '*fastq.gz' | grep -v "Undeter")

sample_dirs=$(for file in $fastq_files; do echo $(basename "$file") | awk '{n=split($0,a,"_");}{sample_dir=a[1]; for(i=2;i<(n-3);i=i+1) {sample_dir= sample_dir"_"a[i];}}END{print sample_dir}'; done | sort | uniq)

for sample_dir in $sample_dirs; do
    mkdir \-p "${output_directory}/${sample_dir}_fastq"
    search_results=$(find "${output_directory}"/"${curr_pipestance_dir}" \-iname *"${sample_dir}"*"fastq.gz")
    echo $search_results
    find "${output_directory}"/"${curr_pipestance_dir}" \-iname *"${sample_dir}"*"fastq.gz" | xargs \-I {} mv {} "${output_directory}"/"${sample_dir}_fastq/"
done
