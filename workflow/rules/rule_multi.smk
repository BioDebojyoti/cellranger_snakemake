# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd
from scripts import create_multi_csv as cmc
from scripts import concat_for_multi as cfm
import shutil
import re

multi_outdir = config_multi["output_multi"]


def generate_multi_input(wc):

    if isinstance(config_multi.get("multi_config_csv"), list):

        config_files_for_multi = config_multi["multi_config_csv"]
        donor_list = []

        for curr_multi_file in config_files_for_multi:

            curr_donor = os.path.basename(curr_multi_file)
            curr_donor = re.sub("_multi_samplesheet.csv", "", curr_donor)
            donor_list.append(curr_donor)

            target_file = f"{curr_donor}_multi_samplesheet.csv"
            destination_path = os.path.join(multi_outdir, target_file)
            # print(target_file)
            # print(destination_path)

            if os.path.abspath(curr_multi_file) != os.path.abspath(destination_path):
                shutil.copyfile(
                    curr_multi_file,
                    destination_path,
                )
    else:

        if isinstance(config_multi.get("input_csvs"), list):
            input_csvs = config_multi["input_csvs"]
        else:
            input_csvs = rules.demultiplex_all.input.fastq_paths

        additional_info_aggr = config_multi["additional_info_aggr"]

        donor_list = cfm.create_donor_dict(
            additional_info_aggr, input_csvs, multi_outdir
        )

        gene_expression_args = (
            {
                "reference": config_multi["gene_expression"]["reference"],
                "probe_set": config_multi["gene_expression"]["probe_set"],
                "filter_probes": config_multi["gene_expression"]["filter_probes"],
                "r1_length": config_multi["gene_expression"]["r1_length"],
                "r2_length": config_multi["gene_expression"]["r2_length"],
                "chemistry": config_multi["gene_expression"]["chemistry"],
                "expect_cells": config_multi["gene_expression"]["expect_cells"],
                "force_cells": config_multi["gene_expression"]["force_cells"],
                "no_secondary": config_multi["gene_expression"]["no_secondary"],
                "create_bam": config_multi["gene_expression"]["create_bam"],
                "check_compatibility": config_multi["gene_expression"][
                    "check_compatibility"
                ],
                "include_introns": config_multi["gene_expression"]["include_introns"],
            }
            if config_multi["gene_expression"]["reference"]
            else None
        )

        vdj_args = (
            {
                "reference": config_multi["vdj"]["reference"],
                "primers": config_multi["vdj"]["primers"],
                "r1_length": config_multi["vdj"]["r1_length"],
                "r2_length": config_multi["vdj"]["r2_length"],
            }
            if config_multi["vdj"]["reference"]
            else None
        )

        for donor_id in donor_list:
            input_multi_csv = f"{multi_outdir}/{donor_id}_multi.csv"
            output_multi_csv = f"{multi_outdir}/{donor_id}_multi_samplesheet.csv"
            # print(output_multi_csv)
            # print(input_multi_csv)

            cmc.create_multi_csv(
                input_csvs=input_multi_csv,
                gene_expression_args=gene_expression_args,
                vdj_args=vdj_args,
                output_file=output_multi_csv,
            )
    return donor_list


donor_list = generate_multi_input("")

# rule cellranger_multi_input:
#     input:
#         csv_files = lambda wc: generate_multi_input


rule cellranger_multi_b4aggr:
    output:
        aggr_input_csv=expand(
            os.path.join("{multi_outdir}", "aggregation_multi.csv"),
            multi_outdir=multi_outdir,
        ),
    input:
        multi_web_summary=expand(
            os.path.join(
                "{multi_outdir}",
                "outs",
                "per_sample_outs",
                "{donor_id}",
                "web_summary.html",
            ),
            donor_id=donor_list,
            multi_outdir=multi_outdir,
        ),
    container:
        "docker://litd/docker-cellranger:v8.0.1"
    params:
        multidir=lambda wc, output: os.path.dirname(output.aggr_input_csv[0]),
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
        csv=os.path.join("{multi_outdir}", "{donor_id}_multi_samplesheet.csv"),
    output:
        multi_web_summary=os.path.join(
            "{multi_outdir}",
            "outs",
            "per_sample_outs",
            "{donor_id}",
            "web_summary.html",
        ),
        multi_count_sample_info=os.path.join(
            "{multi_outdir}",
            "outs",
            "per_sample_outs",
            "{donor_id}",
            "count",
            "sample_molecule_info.h5",
        ),
        multi_count_dir=directory(
            os.path.join(
                "{multi_outdir}",
                "outs",
                "per_sample_outs",
                "{donor_id}",
                "count",
                "sample_filtered_feature_bc_matrix",
            )
        ),
    params:
        id2use=lambda wc: "{donor_id}",
        multi_outdir2use=lambda wc: "{multi_outdir}",
    resources:
        cores=config_count["resources"]["localcores"],
        memory=config_count["resources"]["localmem"],
    container:
        "docker://litd/docker-cellranger:v8.0.1"
    log:
        file=os.path.join("{multi_outdir}", "logs", "multi_{donor_id}.log"),
    benchmark:
        os.path.join("{multi_outdir}", "benchmarks", "benchmarks_{donor_id}_multi.csv")
    shell:
        """   
        cellranger multi \
        --id={params.id2use} \
        --csv={input.csv} \
        --localcores={resources.cores} \
        --localmem={resources.memory} \
        2>&1 | tee -a "{log.file}"; \
        bash scripts/move_pipestance_multi_dir.sh {log.file} {params.multi_outdir2use}; 
        """
