import yaml
import os, sys, argparse
from scripts import create_multi_csv
import subprocess


def get_snakemake_version():
    version_str = subprocess.run(
        ["snakemake", "--version"], capture_output=True, text=True
    ).stdout.strip()
    return tuple(map(int, version_str.split(".")))


snakemake_version = get_snakemake_version()


def get_available_rules(snakemake_version):
    if snakemake_version < (8, 0, 0):
        rules_available = list(rules.__dict__.keys())
    else:
        rules_available = list(rules._rules.keys())
    return rules_available


# Define the function to determine output for cellranger_aggr
# def output_dir_for_cellranger_aggr(wc):

#     rules_available = get_available_rules(snakemake_version)

#     if config_aggr.get("aggr_outdir", None):
#         return config_aggr["aggr_outdir"]
#     elif "cellranger_count" in rules_available:
#         return os.path.join(config_count["output_count"], "aggr_results")
#     elif "cellranger_vdj" in rules_available:
#         return os.path.join(config_vdj["output_vdj"], "aggr_results")
#     else:
#         return os.path.join(config_multi["output_multi"], "aggr_results")


# Define the function to determine input for cellranger_aggr


def input_for_cellranger_aggr(wc):

    rules_available = get_available_rules(snakemake_version)

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
