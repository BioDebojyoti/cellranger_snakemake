# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd


include: os.path.join("rules", "rule_common_config.smk")
include: os.path.join("rules", "rule_mkfastq.smk")
include: os.path.join("rules", "rule_count.smk")
include: os.path.join("rules", "rule_common_aggr.smk")
include: os.path.join("rules", "rule_aggregate.smk")
include: os.path.join("rules", "rule_common_seurat.smk")
include: os.path.join("rules", "rule_seurat.smk")


rule de_results:
    input:
        os.path.join(seurat_outdir, final_rds),
