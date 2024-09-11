import yaml
import sys, argparse
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


# Define the function to determine input for cellranger_count
# def input_for_cellranger_count(wc):
#     # Check if 'fastq path csv or IEM samplesheet' is specified in the config_count file
#     if config_count.get("paths2fastq_file", None):
#         return config_count["paths2fastq_file"]
#     else:
#         # lib_type = rules.cellranger_mkfastq.params.lib_type(wc)
#         # if lib_type in [
#         #     "Gene Expression",
#         #     "gene expression",
#         #     "gene-expression",
#         #     "Gene-Expression",
#         # ]:
#         # If not specified, return the output of rule cellranger_mkfastq
#         return rules.cellranger_mkfastq.output.flag
# else:
#     return None


# Define the function to determine input for cellranger_vdj
# def input_for_cellranger_vdj(wc):
#     # Check if 'fastq path csv or IEM samplesheet' is specified in the config_vdj file
#     if config_vdj.get("paths2fastq4vdj_file", None):
#         return config_vdj["paths2fastq4vdj_file"]
#     else:
#         # If not specified, return the output of rule cellranger_mkfastq
#         return rules.cellranger_mkfastq.output.flag


# # Define the function to determine input for cellranger_aggr
# def input_for_cellranger_aggr(wc):
#     rules_available = list(rules.__dict__.keys())
#     if config_aggr.get("aggr_input_file", None):
#         return config_aggr["aggr_input_file"]
#     elif "cellranger_count" in rules_available:
#         return rules.cellranger_count_b4aggr.output.aggr_input_csv
#     elif "cellranger_vdj" in rules_available:
#         return rules.cellranger_vdj_b4aggr.output.aggr_input_csv
#     else:
#         return rules.cellranger_multi_b4aggr.output.aggr_input_csv


# Define the function to determine output for cellranger_aggr
# def output_dir_for_cellranger_aggr(wc):
#     rules_available = list(rules.__dict__.keys())
#     if config_aggr.get("aggr_outdir", None):
#         return config_aggr["aggr_outdir"]
#     elif "cellranger_count" in rules_available:
#         return config_count["output_count"] + "/aggr_results"
#     elif "cellranger_vdj" in rules_available:
#         return config_vdj["output_vdj"] + "/aggr_results"
#     else:
#         return config_multi["output_multi"] + "/aggr_results"


# # write a multi csv input file
# def get_input_multi_csv():
#     rules_available = list(rules.__dict__.keys())
#     if config_multi.get("multi_csv", None):
#         return config_multi["multi_csv"]
#     else:
#         input_csvs = list(rules.demultiplex_all.input.fastq_paths)
#         gene_expression_args = (
#             {
#                 "reference": config_multi["gene_expression"]["reference"],
#                 "probe_set": config_multi["gene_expression"]["probe_set"],
#                 "filter_probes": config_multi["gene_expression"]["filter_probes"],
#                 "r1_length": config_multi["gene_expression"]["r1_length"],
#                 "r2_length": config_multi["gene_expression"]["r2_length"],
#                 "chemistry": config_multi["gene_expression"]["chemistry"],
#                 "expect_cells": config_multi["gene_expression"]["expect_cells"],
#                 "force_cells": config_multi["gene_expression"]["force_cells"],
#                 "no_secondary": config_multi["gene_expression"]["no_secondary"],
#                 "create_bam": config_multi["gene_expression"]["create_bam"],
#                 "check_compatibility": config_multi["gene_expression"][
#                     "check_compatibility"
#                 ],
#                 "include_introns": config_multi["gene_expression"]["include_introns"],
#             }
#             if config_multi["gene_expression"]["reference"]
#             else None
#         )

#         vdj_args = (
#             {
#                 "reference": config_multi["vdj"]["reference"],
#                 "primers": config_multi["vdj"]["primers"],
#                 "r1_length": config_multi["vdj"]["r1_length"],
#                 "r2_length": config_multi["vdj"]["r2_length"],
#             }
#             if config_multi["vdj"]["reference"]
#             else None
#         )
#         create_multi_csv(
#             input_csvs=input_csvs,
#             gene_expression_args=gene_expression_args,
#             vdj_args=vdj_args,
#             output_file=config_multi["output_multi"]
#             + "/"
#             + config_multi["id"]
#             + "_multi_samplesheet.csv",
#         )
# # Define the function to determine input for cellranger_aggr
# def input_for_cellranger_multi(wc):
#     if config_multi.get("csv", None):
#         return config_aggr["csv"]
#     else:
#         return multi_csv(rules)
