#!/bin/bash

snakemake --profile profile -np --dag | tee workflow.dag.dot
neato -Goverlap=false -Tpng workflow.dag.dot -o dag_cleaned.png
