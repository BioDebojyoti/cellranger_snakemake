# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd


def input_gex_for_seurat(wc):

    rules_available = list(rules.__dict__.keys())

    if config_seurat.get("count_file") is not None:
        return config_seurat["count_file"]
    elif "cellranger_aggr" in rules_available:
        return os.path.join(
            rules.cellranger_aggr.params.aggr_outdir,
            "outs",
            "count",
            "filtered_feature_bc_matrix.h5",
        )
    else:
        return print("Error: File not found or unsupported file type.")


def seurat_input_aggr_csv(wc):

    rules_available = list(rules.__dict__.keys())

    if config_seurat.get("aggr_csv_file") is not None:
        return config_seurat["aggr_csv_file"]
    elif c("cellranger_aggr") in rules_available:
        return rules.cellranger_aggr.input.csv
    else:
        return print("Error: aggregate CSV file missing!!")


def list1_isin_list2(list1, list2):
    present_in_list2 = [elem in list2 for elem in list1]
    if all(present_in_list2):
        return True
    else:
        return False


def vdj_t_flag(wc):
    rules_available = list(rules.__dict__.keys())
    if list1_isin_list2(["cellranger_multi", "cellranger_aggr"], rules_available):
        vdj_t_path = os.path.join(
            rules.cellranger_aggr.params.aggr_outdir,
            "outs",
            "vdj_t",
            "filtered_contig_annotations.csv",
        )
        if os.path.exists(vdj_t_path):
            return f"--vdj-t {vdj_t_path}"
    return ""


def vdj_b_flag(wc):
    # Check if specific rules are available
    rules_available = list(rules.__dict__.keys())
    if list1_isin_list2(["cellranger_multi", "cellranger_aggr"], rules_available):
        vdj_b_path = os.path.join(
            rules.cellranger_aggr.params.aggr_outdir,
            "outs",
            "vdj_b",
            "filtered_contig_annotations.csv",
        )
        if os.path.exists(vdj_b_path):
            return f"--vdj-b {vdj_b_path}"

    return ""


def extra_args_for_seurat(wc, config_seurat):

    # Map config keys to command-line arguments
    arg_mapping = {
        "out_directory": "--out-directory",
        "min_cells": "--min-cells",
        "min_features": "--min-features",
        "max_features": "--max-features",
        "percent_mt": "--percent-mt",
        "percent_rb": "--percent-rb",
        "project": "--project",
        "vdj_t": "--vdj-t",
        "vdj_b": "--vdj-b",
        "layer_column": "--layer-column",
        "condition_column": "--condition-column",
        "integration_method": "--integration-method",
        "enable_SCTransform": "--enable-SCTransform",
        "perform_DE": "--perform-DE",
        "species": "--species",
    }

    # Construct command-line arguments, ensuring missing keys use default values
    additional_args_for_seurat = []
    for param, flag_to_use in arg_mapping.items():
        value = config_seurat.get(param)
        if value is not None:
            # Format logicals as TRUE/FALSE, non-numeric values with quotes, and missing values as NULL
            if isinstance(value, bool):
                formatted_value = "TRUE" if value else "FALSE"
            elif isinstance(value, str):
                formatted_value = f'"{value}"'
            else:
                formatted_value = value

            # Add the formatted argument to the list
            additional_args_for_seurat.append(f"{arg_mapping[param]}={formatted_value}")

    # Join all arguments into a single string
    return " ".join(additional_args_for_seurat)


project_name = config_seurat.get("project") or "singleCell"
final_rds = project_name + "_analysed_seurat.rds"

two_levels_up = os.path.dirname(os.path.dirname(input_gex_for_seurat(None)))
seurat_outdir = config_seurat.get("out_directory") or os.path.join(
    two_levels_up, "seurat_out"
)


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
        extra_args=lambda wc: extra_args_for_seurat(wc, config_seurat),
        vdj_t_args=lambda wc: vdj_t_flag(wc),
        vdj_b_args=lambda wc: vdj_b_flag(wc),
    log:
        file=os.path.join(seurat_outdir, "logs", "seurat.log"),
    shell:
        """
        Rscript scripts/scQCAD.R -f {input.gex} -a {input.aggr_csv} \
        {params.extra_args} {params.vdj_t_args} {params.vdj_b_args} \
        --num-cores={resources.cores} --num-threads={resources.threads} \
        --memory-usage={resources.memory} 2>&1 | tee -a {log.file}; 
        """
