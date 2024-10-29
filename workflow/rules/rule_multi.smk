# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd
from scripts import create_multi_csv as cmc
from scripts import concat_for_multi as cfm
import shutil
import re

multi_outdir = config_multi["output_multi"]

if config_multi.get("multi_config_csv") is not None:
    config_file_for_multi = config_multi["multi_config_csv"]

    print(config_file_for_multi)

    donor_list = os.path.dirname(config_file_for_multi)
    donor_list = donor_list.split("/")[-1]
    donor_list = [re.sub("-", "_", donor_list)]

    target_file = donor_list[0] + "_multi_samplesheet.csv"

    shutil.copyfile(
        config_file_for_multi,
        os.path.join(multi_outdir, target_file),
    )
else:
    input_csvs = config_multi["input_csvs"]

    additional_info_aggr = config_multi["additional_info_aggr"]

    donor_list = cfm.create_donor_dict(additional_info_aggr, input_csvs, multi_outdir)

    # print(donor_list)

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
        input_multi_csv = multi_outdir + "/" + donor_id + "_multi.csv"
        output_multi_csv = multi_outdir + "/" + donor_id + "_multi_samplesheet.csv"
        # print(output_multi_csv)
        # print(input_multi_csv)

        cmc.create_multi_csv(
            input_csvs=input_multi_csv,
            gene_expression_args=gene_expression_args,
            vdj_args=vdj_args,
            output_file=output_multi_csv,
        )


# rule cellranger_multi_input:
#     input:
#         csv_files=rules.demultiplex_all.fastq_paths,


rule cellranger_multi_b4aggr:
    output:
        aggr_input_csv=expand(
            "{multi_outdir}/aggregation_multi.csv", multi_outdir=multi_outdir
        ),
    input:
        multi_web_summary=expand(
            "{multi_outdir}/outs/per_sample_outs/{donor_id}/web_summary.html",
            donor_id=donor_list,
            multi_outdir=multi_outdir,
        ),
    params:
        multidir=multi_outdir,
        additional_info_aggr=config_multi["additional_info_aggr"],
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
        multi_web_summary="{multi_outdir}/outs/per_sample_outs/{donor_id}/web_summary.html",
        # multi_directory_output=directory(
        #     "{multi_outdir}/outs/per_sample_outs/{donor_id}"
        # ),
        multi_count_sample_info="{multi_outdir}/outs/per_sample_outs/{donor_id}/count/sample_molecule_info.h5",
        multi_count_dir=directory(
            "{multi_outdir}/outs/per_sample_outs/{donor_id}/count/sample_filtered_feature_bc_matrix"
        ),
        # multi_vdj_b="{multi_outdir}/outs/per_sample_outs/{donor_id}/vdj_b/vdj_contig_info.pb",
        # multi_vdj_t="{multi_outdir}/outs/per_sample_outs/{donor_id}/vdj_t/vdj_contig_info.pb",
    params:
        id2use="{donor_id}",
        multi_outdir2use="{multi_outdir}",
    resources:
        cores=config_count["resources"]["localcores"],
        memory=config_count["resources"]["localmem"],
    container:
        "docker://litd/docker-cellranger:v8.0.1"
    log:
        file="{multi_outdir}/logs/multi_{donor_id}.log",
    benchmark:
        "{multi_outdir}/benchmarks/benchmarks_{donor_id}_multi.csv"
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
