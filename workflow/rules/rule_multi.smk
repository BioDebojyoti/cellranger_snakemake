# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd


# Configuration
# configfile: os.path.join("config", "config_count.yaml")

multi_outdir = config_count["output_multi"]


# Rule to count features for single library
rule cellranger_multi:
    input:
        id2use=config_multi["id"],
        csv=input_for_cellranger_multi,
    output:
        count="{multi_outdir}/outs/per_sample_outs/{input.id2use}/count/sample_filtered_feature_bc_matrix.h5",
        count_dir=directory(
            "{multi_outdir}/outs/per_sample_outs/{input.id2use}/count/sample_filtered_feature_bc_matrix"
        ),
        vdj_b="{multi_outdir}/outs/per_sample_outs/{input.id2use}/vdj_b/vdj_contig_info.pb",
        vdj_t="{multi_outdir}/outs/per_sample_outs/{input.id2use}/vdj_t/vdj_contig_info.pb",
    params:
        multi_outdir2use="{multi_outdir}",
        add_arguments=config_count["additional_arguments"],
    resources:
        cores=config_count["resources"]["localcores"],
        memory=config_count["resources"]["localmem"],
    container:
        "docker://litd/docker-cellranger:v8.0.1"
    log:
        file="{multi_outdir}/logs/count_{sample}.log",
    benchmark:
        "{multi_outdir}/benchmarks/benchmarks_{input.id2use}_multi.csv"
    shell:
        """   
        cellranger multi \
        --id={input.id2use} \
        --localcores={resources.cores} \
        --localmem={resources.memory} \
        {params.add_arguments} \
        2>&1 | tee -a {log.file}; \
        bash scripts/move_pipestance_multi_dir.sh {log.file} {params.multi_outdir2use}; 
        """
