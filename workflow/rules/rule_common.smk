import yaml
import os, sys, argparse
from scripts import create_multi_csv

# Load config for cellranger_aggr
with open("config/config_aggregate.yaml") as f:
    config_aggr = yaml.safe_load(f)

# Load config for cellranger_mkfastq
with open("config/config_mkfastq.yaml") as f:
    config_mkfastq = yaml.safe_load(f)

# Load config for cellranger_count
with open("config/config_count.yaml") as f:
    config_count = yaml.safe_load(f)

# Load config for cellranger_vdj
with open("config/config_vdj.yaml") as f:
    config_vdj = yaml.safe_load(f)

# Load config for cellranger_multi
with open("config/config_multi.yaml") as f:
    config_multi = yaml.safe_load(f)

# Load config for cellranger_multi
with open("config/config_seurat.yaml") as f:
    config_seurat = yaml.safe_load(f)


# Define the function to determine output for cellranger_aggr
def output_dir_for_cellranger_aggr(wc):
    rules_available = list(rules.__dict__.keys())
    if config_aggr.get("aggr_outdir", None):
        return config_aggr["aggr_outdir"]
    elif "cellranger_count" in rules_available:
        return os.path.join(config_count["output_count"], "aggr_results")
    elif "cellranger_vdj" in rules_available:
        return os.path.join(config_vdj["output_vdj"], "aggr_results")
    else:
        return os.path.join(config_multi["output_multi"], "aggr_results")


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


# Define the function to determine input for cellranger_aggr
def input_for_cellranger_aggr(wc):
    rules_available = list(rules.__dict__.keys())
    if config_aggr.get("aggr_input_file") is not None:
        return config_aggr["aggr_input_file"]
    elif "cellranger_count" in rules_available:
        return rules.cellranger_count_b4aggr.output.aggr_input_csv
    elif "cellranger_vdj" in rules_available:
        return rules.cellranger_vdj_b4aggr.output.aggr_input_csv
    else:
        return rules.cellranger_multi_b4aggr.output.aggr_input_csv


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
