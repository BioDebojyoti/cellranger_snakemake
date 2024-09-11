# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd
from scripts import create_multi_csv
from scripts import create_donor_dict


# Configuration
# configfile: os.path.join("config", "config_count.yaml")

multi_outdir = config_multi["output_multi"]

if config_multi.get("input_csvs", None):
    input_csvs = config_multi["input_csvs"]
else:
    input_csvs = list(rules.demultiplex_all.input.fastq_paths)

additional_info_aggr = config_multi["additional_info_aggr"]

donor_list = create_donor_dict(additional_info_aggr, input_csvs, multi_outdir)


# write a multi csv input file
def get_input_multi_csv(donor_id):

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
    new_multi_csv = (
        config_multi["output_multi"] + "/" + donor_id + "_multi_samplesheet.csv"
    )

    create_multi_csv(
        input_csvs=input_csvs,
        gene_expression_args=gene_expression_args,
        vdj_args=vdj_args,
        output_file=new_multi_csv,
    )

    return new_multi_csv


# Define the function to determine input for cellranger_aggr
# def input_for_cellranger_multi(wc):
#     if config_multi.get("csv", None):
#         return config_aggr["csv"]
#     else:
#         return multi_csv(rules)


# rule cellranger_multi_input:
#     input:
#         csv_files=rules.demultiplex_all.fastq_paths,


rule cellranger_multi_b4aggr:
    output:
        aggr_input_csv=expand(
            "{multi_outdir}/aggregation_multi.csv", multi_outdir=multi_outdir
        ),
    input:
        multi_info_folder=expand(
            directory("{multi_outdir}/outs/per_sample_outs/{donor_id}"),
            donor_id=donor_list,
            multi_outdir=multi_outdir,
        ),
    params:
        multidir=multi_outdir,
        additional_info_aggr=config_count["add_info_aggr"],
    shell:
        """
        bash scripts/get_count_aggr_csv.sh {params.multidir} > {params.multidir}/aggregation_count.csv;
        python scripts/add_info_aggr.py {params.multidir}/aggregation_count.csv {params.additional_info_aggr}
        """


# Rule to count features for single library
rule cellranger_multi:
    input:
        id2use="{donor_id}",
        csv=get_input_multi_csv,
    output:
        multi_count_sample_info="{multi_outdir}/outs/per_sample_outs/{input.id2use}/count/sample_molecule_info.h5",
        multi_count_dir=directory(
            "{multi_outdir}/outs/per_sample_outs/{input.id2use}/count/sample_filtered_feature_bc_matrix"
        ),
        multi_vdj_b="{multi_outdir}/outs/per_sample_outs/{input.id2use}/vdj_b/vdj_contig_info.pb",
        multi_vdj_t="{multi_outdir}/outs/per_sample_outs/{input.id2use}/vdj_t/vdj_contig_info.pb",
    params:
        multi_outdir2use="{multi_outdir}",
    resources:
        cores=config_count["resources"]["localcores"],
        memory=config_count["resources"]["localmem"],
    container:
        "docker://litd/docker-cellranger:v8.0.1"
    log:
        file="{multi_outdir}/logs/multi_{input.id2use}.log",
    benchmark:
        "{multi_outdir}/benchmarks/benchmarks_{input.id2use}_multi.csv"
    shell:
        """   
        cellranger multi \
        --id={input.id2use} \
        --csv={input.csv} \
        --localcores={resources.cores} \
        --localmem={resources.memory} \
        2>&1 | tee -a {log.file}; \
        bash scripts/move_pipestance_multi_dir.sh {log.file} {params.multi_outdir2use}; 
        """
