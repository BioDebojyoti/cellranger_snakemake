# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd


include: os.path.join("rules", "rule_common_config.smk")
include: os.path.join("rules", "rule_mkfastq_count.smk")
include: os.path.join("rules", "rule_count.smk")
include: os.path.join("rules", "rule_common_aggr.smk")
include: os.path.join("rules", "rule_aggregate.smk")
include: os.path.join("rules", "rule_common_seurat.smk")
include: os.path.join("rules", "rule_seurat.smk")


# ruleorder: cellranger_mkfastq > demultiplex_all > cellranger_count > cellranger_count_b4aggr > cellranger_aggr > seurat

# report_file = os.path.join(results_directory, "report.html")


rule de_results:
    input:
        os.path.join(seurat_outdir, final_rds),


onsuccess:
    print("Workflow completed successfully!")
    shell(f"snakemake --report {report_file}")
