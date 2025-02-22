# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd
from scripts import create_multi_csv as cmc
from scripts import concat_for_multi as cfm
import shutil
import re


report: "../report/workflow.rst"


multi_outdir = os.path.join(results_directory, "multi_results")

donor_list = list(
    set(pd.read_csv(config_multi["additional_info_aggr"])["donor"].tolist())
)


checkpoint cellranger_multi_input_prep:
    output:
        expand(
            os.path.join(multi_outdir, "{donor_id}_multi_samplesheet.csv"),
            donor_id=donor_list,
        ),
    input:
        fastq_paths=lambda wc: list(set(rules.demultiplex_all.input.csvs)),
    params:
        config_file="config/config_multi.yaml",
        multi_outdir=lambda wc, output: os.path.dirname(output[0]),
    log:
        os.path.join(multi_outdir, "logs", "multi_input_prep.log"),
    container:
        "docker://biodebojyoti/pythoncellranger:v1.0.0"
    shell:
        """
        python scripts/get_multi_input.py {params.config_file} {params.multi_outdir} {input.fastq_paths}
        """


rule cellranger_multi_b4aggr:
    output:
        aggr_input_csv=os.path.join(multi_outdir, "aggregation_multi.csv"),
    input:
        multi_count_sample_info=expand(
            os.path.join(
                multi_outdir,
                "{donor_id}",
                "outs",
                "per_sample_outs",
                "{donor_id}",
                "count",
                "sample_molecule_info.h5",
            ),
            donor_id=donor_list,
        ),
    container:
        "docker://biodebojyoti/pythoncellranger:v1.0.0"
    # conda:
    #     "../envs/minimal_python.yaml"
    params:
        multidir=lambda wc: multi_outdir,
        additional_info_aggr=config_multi["additional_info_aggr"],
    log:
        os.path.join(multi_outdir, "logs", "multi_pre_aggr.log"),
    shell:
        """
        bash scripts/get_multi_aggr_csv.sh {params.multidir} > {output.aggr_input_csv};
        python scripts/add_info_aggr_multi.py {output.aggr_input_csv} {params.additional_info_aggr}
        """


# Rule to count features for single library
rule cellranger_multi:
    input:
        csv=os.path.join(multi_outdir, "{donor_id}_multi_samplesheet.csv"),
    output:
        multi_web_summary=report(
            directory(
                os.path.join(
                    multi_outdir, "{donor_id}", "outs", "per_sample_outs", "{donor_id}"
                )
            ),
            htmlindex="web_summary.html",
            caption="../report/cellranger_multi_mapping.rst",
            category="cellranger_multi",
            subcategory="{donor_id}",
            labels={"donor": "{donor_id}"},
        ),
        multi_count_sample_info=os.path.join(
            multi_outdir,
            "{donor_id}",
            "outs",
            "per_sample_outs",
            "{donor_id}",
            "count",
            "sample_molecule_info.h5",
        ),
        multi_count_dir=directory(
            os.path.join(
                multi_outdir,
                "{donor_id}",
                "outs",
                "per_sample_outs",
                "{donor_id}",
                "count",
                "sample_filtered_feature_bc_matrix",
            )
        ),
    params:
        id2use="{donor_id}",
        multi_outdir2use=lambda wc: os.path.join(multi_outdir, f"{wc.donor_id}"),
    resources:
        cores=config_count["resources"]["localcores"],
        memory=config_count["resources"]["localmem"],
    container:
        "docker://litd/docker-cellranger:v8.0.1"
    log:
        file=os.path.join(multi_outdir, "logs", "multi_{donor_id}.log"),
    benchmark:
        os.path.join(multi_outdir, "benchmarks", "benchmarks_{donor_id}_multi.csv")
    shell:
        """   
        cellranger multi \
        --id={params.id2use} \
        --csv={input.csv} \
        --localcores={resources.cores} \
        --localmem={resources.memory} \
        >> {log.file} 2>&1; \
        bash scripts/move_pipestance_multi_dir.sh {log.file} {params.multi_outdir2use}; 
        """
