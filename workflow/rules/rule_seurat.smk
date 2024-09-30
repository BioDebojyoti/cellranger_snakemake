# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd


# Configuration
configfile: os.path.join("config", "config_seurat.yaml")


def input_gex_for_seurat(wc):

    rules_available = list(rules.__dict__.keys())

    if config_seurat.get("count_matrix_file", None):
        return config_seurat["count_matrix_file"]
    if c("cellranger_multi", "cellranger_aggr") in rules_available:
        return os.path.join(
            rules.cellranger_aggr.params.aggr_outdir,
            "outs",
            "count",
            "filtered_feature_bc_matrix.h5",
        )

def vdj_t_flag(wc):
    
    rules_available = list(rules.__dict__.keys())
    if c("cellranger_multi", "cellranger_aggr") in rules_available:
        if os.path.exists(
            rules.cellranger_aggr.params.aggr_outdir
            + "/outs/vdj_t/filtered_contig_annotations.csv"
        ):
            return "--vdj-t " os.path.join(
                    rules.cellranger_aggr.params.aggr_outdir,
                    "outs",
                    "vdj_t",
                    "filtered_contig_annotations.csv",
                ) 
        else:
        return ""

def vdj_b_flag(wc):
    
    rules_available = list(rules.__dict__.keys())
    if c("cellranger_multi", "cellranger_aggr") in rules_available:
        if os.path.exists(
            rules.cellranger_aggr.params.aggr_outdir
            + "/outs/vdj_b/filtered_contig_annotations.csv"
        ):
            return "--vdj-b " os.path.join(
                    rules.cellranger_aggr.params.aggr_outdir,
                    "outs",
                    "vdj_b",
                    "filtered_contig_annotations.csv",
                ) 
        else:
        return ""


# Rule to aggregate libraries (optional)
rule seurat:
    input:
        gex =le,
    output:
        seurat_rds="{seurat_outdir}/",
    resources:
        cores=config["resources"]["localcores"],
        memory=config["resources"]["localmem"],
    container:
        "docker://biodebojyoti/crossplatformcellranger:8.0.1"
    params:
        aggr_id=config["aggregation_id"],
        aggr_outdir="{count_outdir}/aggr_results",
    log:
        file="{count_outdir}/logs/aggr.log",
    shell:
        """
        cellranger aggr --id={params.aggr_id} --csv={input.csv} --normalize=mapped \
        --localcores={resources.cores} \
        --localmem={resources.memory} \
        2>&1 | tee -a {log.file};
        bash scripts/move_pipestance_count_dir.sh {log.file} {params.aggr_outdir}; 
        """
