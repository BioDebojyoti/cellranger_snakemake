#!/bin/bash

dag_type=$2
workflow_dot_file=$1
workflow_dot_file=$dag_type"_"$workflow_dot_file

snakemake --profile profile -np --$dag_type | tee $workflow_dot_file
awk 'BEGIN{print "@startuml"}{print $0}END{print "@enduml"}' $workflow_dot_file >${workflow_dot_file%.dot}".puml"
