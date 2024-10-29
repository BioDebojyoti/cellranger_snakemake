# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd

# Configuration
# configfile: os.path.join("config", "config_aggregate.yaml")

aggr_outdir = output_dir_for_cellranger_aggr(" ")


# Define the function to determine input for cellranger_aggr
def input_for_cellranger_aggr(wc):
    rules_available = list(rules.__dict__.keys())
    if config_aggr.get("aggr_input_file") is not None:
        return config_aggr["aggr_input_file"]
    elif "cellranger_count" in rules_available:
        return rules.cellranger_count_b4aggr.output.aggr_input_csv
    elif "cellranger_vdj" in rules_available:
        return rules.cellranger_vdj_b4aggr.output.aggr_input_csv
    else:
        return rules.cellranger_multi_b4aggr.output.aggr_input_csv


# Rule to aggregate libraries
rule cellranger_aggr:
    input:
        csv=input_for_cellranger_aggr,
        # aggr_input_file if aggr_input_file != None else rules.b4all.output.aggr_input_csv
    output:
        aggr_h5="{aggr_outdir}/outs/count/filtered_feature_bc_matrix.h5",
    resources:
        cores=config_aggr["resources"]["localcores"],
        memory=config_aggr["resources"]["localmem"],
    container:
        "docker://litd/docker-cellranger:v8.0.1"
    params:
        aggr_id=config_aggr["aggregation_id"],
        aggr_outdir="{aggr_outdir}",
        aggr_normalize=config_aggr.get("normalize") or "none",
    log:
        file="{aggr_outdir}/logs/aggr.log",
    benchmark:
        "{aggr_outdir}/benchmarks/benchmark_{params.aggr_id}.csv"
    shell:
        """
        cellranger aggr --id={params.aggr_id} --csv={input.csv} --normalize={params.aggr_normalize} \
        --localcores={resources.cores} \
        --localmem={resources.memory} \
        2>&1 | tee -a {log.file};
        bash scripts/move_pipestance_aggr_dir.sh {log.file} {params.aggr_outdir}; 
        """
