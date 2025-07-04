# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd


include: os.path.join("rules", "rule_common_config.smk")
include: os.path.join("rules", "rule_mkfastq_multi.smk")
include: os.path.join("rules", "rule_multi.smk")
include: os.path.join("rules", "rule_common_aggr.smk")
include: os.path.join("rules", "rule_aggregate.smk")
include: os.path.join("rules", "rule_common_seurat.smk")
include: os.path.join("rules", "rule_seurat.smk")


ruleorder: cellranger_mkfastq > demultiplex_all > cellranger_multi_input_prep > cellranger_multi > cellranger_multi_b4aggr > cellranger_aggr > seurat_step


rule de_results:
    input:
        os.path.join(seurat_outdir, final_rds),


onsuccess:
    print("Workflow completed successfully!")
    shell(f"snakemake -s Snakefile_multi_full.smk --report {report_file}")
