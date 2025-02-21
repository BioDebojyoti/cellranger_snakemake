# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd

aggr_outdir = os.path.join(results_directory, "aggr_results")


# Rule to aggregate libraries
rule cellranger_aggr:
    input:
        csv=input_for_cellranger_aggr,
    output:
        aggr_h5=os.path.join(
            "{aggr_outdir}", "outs", "count", "filtered_feature_bc_matrix.h5"
        ),
        WebSummary=report(
            directory(os.path.join("{aggr_outdir}", "outs")),
            htmlindex="web_summary.html",
            caption="../report/cellranger_aggr.rst",
            category="cellranger_aggr",
            labels={"ID used": config_aggr["aggregation_id"]},
        ),
    resources:
        cores=config_aggr["resources"]["localcores"],
        memory=config_aggr["resources"]["localmem"],
    container:
        "docker://litd/docker-cellranger:v8.0.1"
    params:
        aggr_id=config_aggr["aggregation_id"],
        aggr_outdir=aggr_outdir,
        aggr_normalize=config_aggr.get("normalize") or "none",
    log:
        file=os.path.join("{aggr_outdir}", "logs", "aggr.log"),
    benchmark:
        os.path.join("{aggr_outdir}", "benchmarks", "benchmark_{params.aggr_id}.csv")
    shell:
        """
        cellranger aggr --id={params.aggr_id} --csv={input.csv} --normalize={params.aggr_normalize} \
        --localcores={resources.cores} \
        --localmem={resources.memory} \
        >> {log.file} 2>&1; \
        bash scripts/move_pipestance_aggr_dir.sh {log.file} {params.aggr_outdir}; 
        """
