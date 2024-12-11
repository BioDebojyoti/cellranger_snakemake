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
sample_to_run = {
    k: v for v in run_bcl_sample_dict.keys() for k in run_bcl_sample_dict[v]
}

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

bcl_run_index_dir_dict = {}

for run, fastq_dir_path in sample_directories_dict_transformed.values():
    if run not in bcl_run_index_dir_dict:
        bcl_run_index_dir_dict[run] = []
    bcl_run_index_dir_dict[run].append(fastq_dir_path)


# print(bcl_run_index_dir_dict)

directories = bcl_run_index_dir_dict.values()

df = pd.DataFrame.from_dict(
    {
        key: {"sample": key, "bcl_run_index": value[0], "fastq_dir_path": value[1]}
        for key, value in sample_directories_dict_transformed.items()
    },
    orient="index",
).reset_index(drop=True)

# Reorder columns to match the desired order
df = df[["sample", "bcl_run_index", "fastq_dir_path"]]
df["fastq_outdirectory"] = [os.path.dirname(d) for d in df["fastq_dir_path"].tolist()]


flag_files = []
for index, row in df.iterrows():
    bcl = row["bcl_run_index"]
    curr_flag_file = os.path.join(
        row["fastq_outdirectory"],
        f"mkfastq_success_{bcl}.csv",
    )
    if curr_flag_file not in flag_files:
        flag_files.append(curr_flag_file)

flag_dictionary = {}
for i, v in enumerate(df["bcl_run_index"].unique()):
    curr_df = df[df["bcl_run_index"] == v].copy()
    curr_flag_file = os.path.join(
        curr_df["fastq_outdirectory"].unique()[0],
        f"mkfastq_success_{bcl}.csv",
    )
    flag_dictionary[v] = [curr_flag_file, curr_df["fastq_dir_path"].tolist()]

flag_files = [f[0] for f in flag_dictionary.values()]
fastq_folders = [directory(d) for f in flag_dictionary.values() for d in f[1]]
fastq_outdirectory = df["fastq_outdirectory"].unique()[0]


rule demultiplex_all:
    input:
        expand(
            os.path.join(fastq_outdirectory, "mkfastq_success_{bcl_run_index}.csv"),
            bcl_run_index=bcl_run_indexes,
        ),
        # expand(os.path.join(fastq_outdirectory, "{sample}_fastq"), sample=samples),


checkpoint fastq_folder_collect:
    output:
        directory(
            os.path.join(fastq_outdirectory, "{sample}_fastq"),
        ),
    input:
        lambda wc: os.path.join(
            fastq_outdirectory, f"mkfastq_success_{sample_to_run[wc.sample]}.csv"
        ),


# input:
#     lambda wc: os.path.join(
#         "{fastq_outdirectory}", "mkfastq_success_{bcl_run_index}.csv"
#     ),


# Rule to generate FASTQ files if starting from BCL files
checkpoint cellranger_mkfastq:
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
        # "cellranger.v8.0.1.sif"
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
        bash scripts/move_pipestance_mkfastq_dir.sh {log:q} {params.outdir2use:q};
        bash scripts/get_fastq_csv.sh {params.outdir2use:q} "{params.lib_type}" > {output.flag};
        """
