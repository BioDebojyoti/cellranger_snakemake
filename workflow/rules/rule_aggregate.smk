# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd

# Configuration
# configfile: os.path.join("config", "config_aggregate.yaml")

aggr_outdir = output_dir_for_cellranger_aggr(" ")


# Rule to aggregate libraries (optional)
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
    log:
        file="{aggr_outdir}/logs/aggr.log",
    shell:
        """
        cellranger aggr --id={params.aggr_id} --csv={input.csv} --normalize=mapped \
        --localcores={resources.cores} \
        --localmem={resources.memory} \
        2>&1 | tee -a {log.file};
        bash scripts/move_pipestance_count_dir.sh {log.file} {params.aggr_outdir}; 
        """
