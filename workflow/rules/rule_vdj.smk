# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd


# Configuration
# configfile: os.path.join("config", "config_vdj.yaml")


# fastq_df_vdj = pd.read_csv(config_vdj["paths2fastq4vdj_file"])
fastq_df_vdj = pd.read_csv(input_for_cellranger_vdj(" "))
vdj_samples_seqs = fastq_df_vdj["sample"].tolist()
vdj_samples_paths = fastq_df_vdj["fastq"].tolist()
vdj_fastq_dict = {s: vdj_samples_paths[i] for i, s in enumerate(vdj_samples_seqs)}


vdj_outdir = config_vdj["output_vdj"]

vdj_cores = config_vdj["resources"]["localcores"]
vdj_memory = config_vdj["resources"]["localmem"]

# Define the upper limits
max_cores = config_vdj["resources"]["max_cores"]
max_memory = config_vdj["resources"]["max_memory"]


rule cellranger_vdj_b4aggr:
    output:
        aggr_input_csv=expand(
            "{vdj_outdir}/aggregation_count.csv", vdj_outdir=vdj_outdir
        ),
    input:
        # dummy inputs to create csv for cellranger_aggr  
        # flag=os.path.join(outdir, "mkfastq.sucess.csv"),
        count_info_h5=expand(
            "{vdj_outdir}/{sample}_count/outs/molecule_info.h5",
            sample=list(vdj_fastq_dict.keys()),
            vdj_outdir=vdj_outdir,
        ),
    params:
        vdjdir=vdj_outdir,
        additional_info_aggr=config_vdj["add_info_aggr"],
    shell:
        """
        bash scripts/get_vdj_aggr_csv.sh {params.vdjdir} > {params.vdjdir}/aggregation_count.csv;
        python scripts/add_info_aggr.py {params.vdjdir}/aggregation_count.csv {params.additional_info_aggr}
        """


# Rule to vdj for single library
rule cellranger_vdj:
    input:
        fastq_folder=lambda wc: vdj_fastq_dict[wc.sample],
        vdj_reference=config_vdj["vdj_reference"],
    output:
        annotations_bed="{vdj_outdir}/{sample}_vdj/outs/all_contig_annotations.bed",
        annotations_json="{vdj_outdir}/{sample}_vdj/outs/all_contig_annotations.json",
        metrics_summary_csv="{vdj_outdir}/{sample}_vdj/outs/metrics_summary.csv",
        vdj_contig_info_pb="{vdj_outdir}/{sample}_vdj/outs/vdj_contig_info.pb",
        vdj_reference=directory("{vdj_outdir}/{sample}_vdj/outs/vdj_reference"),
        vdj_cloupe="{vdj_outdir}/{sample}_vdj/outs/vloupe.vloupe",
        vdj_web_summary="{vdj_outdir}/{sample}_vdj/outs/web_summary.html",
    params:
        id2use="pipestance_{sample}",
        sample="{sample}",
        vdj_outdir2use="{vdj_outdir}/{sample}_vdj",
        add_arguments=config_vdj["additional_arguments"],
    resources:
        cores=lambda wc, attempt: min(vdj_cores * attempt, max_cores),
        memory=lambda wc, attempt: min(vdj_memory * attempt, max_memory),
    container:
        "docker://litd/docker-cellranger:v8.0.1"
    log:
        file="{vdj_outdir}/logs/vdj_{sample}.log",
    benchmark:
        "{vdj_outdir}/benchmarks/benchmarks_{sample}_vdj.csv"
    shell:
        """   
        cellranger vdj \
        --id={params.id2use} \
        --reference={input.vdj_reference} \
        --fastqs={input.fastq_folder} \
        --sample={params.sample} \
        --localcores={resources.cores} \
        --localmem={resources.memory} \
        {params.add_arguments} \
        2>&1 | tee -a {log.file}; \
        bash scripts/move_pipestance_vdj_dir.sh {log.file} {params.vdj_outdir2use}; 
        """
