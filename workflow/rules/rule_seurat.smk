# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd


project_name = config_seurat.get("project") or "singleCell"
final_rds = f"{project_name}_analysed_seurat.rds"

seurat_outdir = os.path.join(fastq_outdirectory, "seurat_out")


# Rule to run seurat on cellranger output
rule seurat:
    input:
        gex=input_gex_for_seurat,
        aggr_csv=seurat_input_aggr_csv,
    output:
        seurat_rds=os.path.join(seurat_outdir, final_rds),
    resources:
        cores=config_seurat["resources"]["cores"],
        threads=config_seurat["resources"]["threads"],
        memory=config_seurat["resources"]["memory"],
    conda:
        "d_seurat510"
    params:
        seurat_out_directory=seurat_outdir,
        extra_args=lambda wc: extra_args_for_seurat(wc, config_seurat),
        vdj_t_args=lambda wc: vdj_t_flag(wc),
        vdj_b_args=lambda wc: vdj_b_flag(wc),
    log:
        file=os.path.join(seurat_outdir, "logs", "seurat.log"),
    shell:
        """
        Rscript scripts/scQCAD.R -f {input.gex} -a {input.aggr_csv} \
        {params.extra_args} {params.vdj_t_args} {params.vdj_b_args} \
        --out-directory={params.seurat_out_directory} \
        --num-cores={resources.cores} --num-threads={resources.threads} \
        --memory-usage={resources.memory} 2>&1 | tee -a {log.file}; 
        """
