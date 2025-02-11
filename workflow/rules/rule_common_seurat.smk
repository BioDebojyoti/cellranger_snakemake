import yaml
import os, sys, argparse
from scripts import create_multi_csv


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
    else:
        if config_aggr.get("aggr_input_file") is not None:
            return config_aggr["aggr_input_file"]
        elif "cellranger_count_b4aggr" in rules_available:
            return rules.cellranger_count_b4aggr.output.aggr_input_csv
        elif "cellranger_multi_b4aggr" in rules_available:
            return rules.cellranger_multi_b4aggr.output.aggr_input_csv
        elif "cellranger_vdj_b4aggr" in rules_available:
            return rules.cellranger_vdj_b4aggr.output.aggr_input_csv
        else:
            print("Something went wrong!!")


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
        # "out_directory": "--out-directory",
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
