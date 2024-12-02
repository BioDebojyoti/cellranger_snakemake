# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd

count_outdir = config_count["output_count"]

donor_list = list(set(pd.read_csv(config_count["add_info_aggr"])["donor"].tolist()))


def samples_for_count():
    rules_available = list(rules.__dict__.keys())
    if "demultiplex_all" in rules_available:
        samples_to_process = samples
        global fastq_dir_paths
        sample_count_fastq_dict = {
            s: p for s, p in zip(samples_to_process, fastq_dir_paths)
        }
    else:
        samples_to_process = pd.read_csv(config_count["paths2fastq_file"])
        samples_to_process = samples_to_process_df["sample"].tolist()
        fastq_dir_paths = samples_to_process_df["fastq"].tolist()
        sample_count_fastq_dict = {
            s: p for s, p in zip(samples_to_process, fastq_dir_paths)
        }
    return [samples_to_process, sample_count_fastq_dict]


dict_elements = samples_for_count()
samples_to_process = dict_elements[0]
sample_count_fastq_dict = dict_elements[1]

# print(samples_to_process)
# print(sample_count_fastq_dict)

# def create_sample_dict_from_csv(file_path):
#     df = pd.read_csv(file_path, index_col=0)
#     sample_dict = df.apply(lambda row: row.dropna().tolist(), axis=1).to_dict()
#     return [sample_dict.keys(), sample_dict.values(), sample_dict]


def load_fastq_dict(fastq_dict_file):

    df = pd.read_csv(fastq_dict_file, index_col=0)
    return df.to_dict(orient="index")


# Dynamically define fastq_file_path based on user input or mkfastq output
# fastq_file_path = (
#     config_count["paths2fastq_file"]
#     if config_count["paths2fastq_file"] is not None
#     else expand(
#         os.path.join("{fastq_outdirectory}", "mkfastq_success_{bcl_run_index}.csv"),
#         zip,
#         fastq_outdirectory=list(fastq_outdirectory_dict.values()),
#         bcl_run_index=bcl_run_index_gex,
#     )
# )


# Rule to create the FASTQ dictionary dynamically
# rule create_fastq_dict:
#     input:
#         samplesheet_or_fastq_files=lambda: config_count["paths2fastq_file"]
#         or rules.demultiplex_all.input.fastq_paths,
#     output:
#         fastq_dict=os.path.join(config_count["output_count"], "fastq_dict.csv"),
#     params:
#         additional_info=config_count["add_info_aggr"],
#     container:
#         "docker://litd/docker-cellranger:v8.0.1"
#     log:
#         os.path.join(count_outdir, "logs", "count_fastq_dict.log"),
#     shell:
#         """
#         python scripts/create_count_dict.py {params.additional.info} \
#         {input.samplesheet_or_fastq_files} > {output.fastq_dict}
#         """


rule cellranger_count_b4aggr:
    output:
        aggr_input_csv=expand(
            os.path.join("{count_outdir}", "aggregation_count.csv"),
            count_outdir=count_outdir,
        ),
    input:
        # flag=os.path.join(outdir, "mkfastq.sucess.csv"),
        count_info_h5=expand(
            os.path.join(
                "{count_outdir}", "{sample}_count", "outs", "molecule_info.h5"
            ),
            sample=samples_to_process,
            count_outdir=count_outdir,
        ),
    log:
        os.path.join(count_outdir, "logs", "count_aggr.log"),
    params:
        countdir=lambda wc, output: os.path.dirname(output.aggr_input_csv[0]),
        additional_info_aggr=config_count["add_info_aggr"],
    container:
        "docker://litd/docker-cellranger:v8.0.1"
    shell:
        """
        bash scripts/get_count_aggr_csv.sh {params.countdir} > {params.countdir}/aggregation_count.csv;
        python scripts/add_info_aggr.py {params.countdir}/aggregation_count.csv {params.additional_info_aggr}
        """


# Rule to count features for single library
rule cellranger_count:
    input:
        # fastq_folder=rules.create_fastq_dict.output,
        fastq_folder=lambda wc: sample_count_fastq_dict[wc.sample],
        transcriptome=config_count["transcriptome"],
    output:
        Permoleculereadinformation=os.path.join(
            "{count_outdir}", "{sample}_count", "outs", "molecule_info.h5"
        ),
        FilteredBCmatricesHDF5=os.path.join(
            "{count_outdir}", "{sample}_count", "outs", "filtered_feature_bc_matrix.h5"
        ),
    params:
        id2use=lambda wc: f"pipestance_{wc.sample}",
        sample=lambda wc: f"{wc.sample}",
        count_outdir2use=lambda wc: f"{count_outdir}/{wc.sample}_count",
        add_arguments=config_count["additional_arguments"],
    resources:
        cores=config_count["resources"]["localcores"],
        memory=config_count["resources"]["localmem"],
    container:
        "docker://litd/docker-cellranger:v8.0.1"
    log:
        file=os.path.join("{count_outdir}", "logs", "count_{sample}.log"),
    benchmark:
        os.path.join("{count_outdir}", "benchmarks", "benchmarks_{sample}_count.csv")
    shell:
        """   
        cellranger count \
        --id={params.id2use} \
        --transcriptome={input.transcriptome} \
        --fastqs={input.fastq_folder} \
        --sample={params.sample} \
        --localcores={resources.cores} \
        --localmem={resources.memory} \
        {params.add_arguments} \
        2>&1 | tee -a {log.file}; \
        bash scripts/move_pipestance_count_dir.sh {log.file} {params.count_outdir2use}; 
        """
