# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd


include: os.path.join("rules", "rule_common.smk")
include: os.path.join("rules", "rule_mkfastq.smk")
include: os.path.join("rules", "rule_count.smk")
include: os.path.join("rules", "rule_aggregate.smk")
include: os.path.join("rules", "rule_seurat.smk")


# print(list(rules.__dict__.keys()))


rule de_results:
    input:
        os.path.join(seurat_outdir, final_rds),
