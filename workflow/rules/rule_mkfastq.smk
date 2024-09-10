# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys, re
import pandas as pd

# Configuration
# configfile: os.path.join("config", "config_mkfastq.yaml")

# Define paths based on configuration
# bcl_folder = config_mkfastq["bcl_folder"]

# outdir = config_mkfastq["outputdir"]
# library_type = config_mkfastq["library_type"]

add_args = [
    (
        config_mkfastq["additional_arguments"]
        if config_mkfastq["additional_arguments"] != None
        else ""
    )
]

mkfastq_cores = config_mkfastq["resources"]["localcores"]
mkfastq_memory = config_mkfastq["resources"]["localmem"]

# Define the upper limits
max_cores = config_mkfastq["resources"]["max_cores"]
max_memory = config_mkfastq["resources"]["max_memory"]


combinations_df = pd.read_csv(config_mkfastq["bcl_folder_paths"])

run_bcl_paths_dict = dict(
    zip(combinations_df["bcl_run_index"], combinations_df["run_bcl_path"])
)

fastq_outdirectory_dict = dict(
    zip(combinations_df["bcl_run_index"], combinations_df["fastq_outdirectory"])
)

feature_type_dict = dict(
    zip(combinations_df["bcl_run_index"], combinations_df["feature_type"])
)

iem_samplesheet_dict = dict(
    zip(combinations_df["bcl_run_index"], combinations_df["iem_samplesheet"])
)

samplesheet_4_bcl_dict = dict(
    zip(combinations_df["bcl_run_index"], combinations_df["samplesheet_4_bcl"])
)


rule demultiplex_all:
    input:
        fastq_paths=expand(
            "{fastq_outdirectory}/mkfastq_success_{bcl_run_index}.csv",
            zip,
            fastq_outdirectory=list(fastq_outdirectory_dict.values()),
            bcl_run_index=list(fastq_outdirectory_dict.keys()),
        ),


# Rule to generate FASTQ files if starting from BCL files
rule cellranger_mkfastq:
    input:
        bcl_path=lambda wc: run_bcl_paths_dict[wc.bcl_run_index],
    output:
        flag="{fastq_outdirectory}/mkfastq_success_{bcl_run_index}.csv",
    resources:
        cores=lambda wc, attempt: min(mkfastq_cores * attempt, max_cores),
        memory=lambda wc, attempt: min(mkfastq_memory * attempt, max_memory),
    params:
        sampleinfo=lambda wc: (
            " --samplesheet=" + samplesheet_4_bcl_dict[wc.bcl_run_index]
            if iem_samplesheet_dict[wc.bcl_run_index]
            else " --csv=" + samplesheet_4_bcl_dict[wc.bcl_run_index]
        ),
        args2add=add_args[0],
        outdir2use=lambda wc: fastq_outdirectory_dict[wc.bcl_run_index],
        lib_type=lambda wc: feature_type_dict[wc.bcl_run_index],
    container:
        "docker://litd/docker-cellranger:v8.0.1"
    log:
        os.path.join("{fastq_outdirectory}/logs/mkfastq_{bcl_run_index}.log"),
    benchmark:
        os.path.join(
            "{fastq_outdirectory}/benchmarks/benchmark_mkfastq_{bcl_run_index}"
        )
    shell:
        """
        cellranger mkfastq --run={input.bcl_path} \
        --output-dir={params.outdir2use} \
        {params.sampleinfo} \
        {params.args2add} \
        --localcores={resources.cores} \
        --localmem={resources.memory} \
        2>&1 | tee -a {log};
        bash scripts/move_pipestance_mkfastq_dir.sh {log} {params.outdir2use};
        bash scripts/get_fastq_csv.sh {params.outdir2use} "{params.lib_type}" > {output.flag};
        """
