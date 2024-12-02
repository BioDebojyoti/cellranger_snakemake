# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys, re
import pandas as pd


add_args = [
    (
        config_mkfastq["additional_arguments"]
        if config_mkfastq["additional_arguments"] is not None
        else ""
    )
]

mkfastq_cores = config_mkfastq["resources"]["localcores"]
mkfastq_memory = config_mkfastq["resources"]["localmem"]

# Define the upper limits
max_cores = config_mkfastq["resources"]["max_cores"]
max_memory = config_mkfastq["resources"]["max_memory"]


combinations_df = pd.read_csv(config_mkfastq["bcl_folder_paths"])


dfs = pd.DataFrame()

for i, run_index in enumerate(combinations_df["bcl_run_index"].tolist()):
    if combinations_df.iloc[i]["iem_samplesheet"]:
        df_curr_run = extract_data_block(combinations_df.iloc[i]["samplesheet_4_bcl"])
        df_curr_run.columns.values[1] = "Sample"
    else:
        df_curr_run = pd.read_csv(combinations_df.iloc[i]["samplesheet_4_bcl"])

    df_curr_run["bcl_run_index"] = combinations_df.iloc[i]["bcl_run_index"]
    # Concatenate the current DataFrame with the overall DataFrame
    dfs = pd.concat([dfs, df_curr_run], ignore_index=True)

run_bcl_sample_dict = dfs.groupby("bcl_run_index")["Sample"].apply(list).to_dict()


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

bcl_run_index_gex = [
    i
    for i, v in feature_type_dict.items()
    if v in ["Gene Expression", "gene expression", "gene-expression", "Gene-Expression"]
]

sample_directories_dict_transformed = {
    sample: (
        bcl_run_index,
        os.path.join(fastq_outdirectory_dict[bcl_run_index], f"{sample}_fastq"),
    )
    for bcl_run_index in run_bcl_sample_dict.keys()
    for sample in run_bcl_sample_dict[bcl_run_index]
}


samples = list(sample_directories_dict_transformed.keys())
bcl_run_indexes = [values[0] for values in sample_directories_dict_transformed.values()]
fastq_dir_paths = [values[1] for values in sample_directories_dict_transformed.values()]

# print(sample_directories_dict_transformed)


rule demultiplex_all:
    input:
        fastq_paths=expand(
            os.path.join(
                "{fastq_outdirectory}", "mkfastq_success_{bcl_run_index}.csv"
            ),
            zip,
            fastq_outdirectory=list(fastq_outdirectory_dict.values()),
            bcl_run_index=list(fastq_outdirectory_dict.keys()),
        ),
    output:
        demultiplexed_samples=expand(
            "{directory}",
            directory=[
                sample_directories_dict_transformed[sample][1]
                for sample in sample_directories_dict_transformed.keys()
            ],
        ),


# Rule to generate FASTQ files if starting from BCL files
rule cellranger_mkfastq:
    input:
        bcl_path=lambda wc: run_bcl_paths_dict[wc.bcl_run_index],
    output:
        flag=os.path.join("{fastq_outdirectory}", "mkfastq_success_{bcl_run_index}.csv"),
    resources:
        cores=lambda wc, attempt: min(mkfastq_cores * attempt, max_cores),
        memory=lambda wc, attempt: min(mkfastq_memory * attempt, max_memory),
    params:
        sampleinfo=lambda wc: (
            f" --samplesheet={samplesheet_4_bcl_dict[wc.bcl_run_index]}"
            if iem_samplesheet_dict[wc.bcl_run_index]
            else f" --csv={samplesheet_4_bcl_dict[wc.bcl_run_index]}"
        ),
        args2add=add_args[0],
        outdir2use=lambda wc: fastq_outdirectory_dict[wc.bcl_run_index],
        lib_type=lambda wc: feature_type_dict[wc.bcl_run_index],
    container:
        "docker://litd/docker-cellranger:v8.0.1"
    log:
        os.path.join("{fastq_outdirectory}", "logs", "mkfastq_{bcl_run_index}.log"),
    benchmark:
        os.path.join(
            "{fastq_outdirectory}", "benchmarks", "benchmark_mkfastq_{bcl_run_index}"
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
