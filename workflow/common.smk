import yaml

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


# Define the function to determine input for cellranger_aggr
def input_for_cellranger_aggr(wc):
    rules_available = list(rules.__dict__.keys())
    if config_aggr.get("aggr_input_file", None):
        return config_aggr["aggr_input_file"]
    elif "cellranger_count" in rules_available:
        return rules.cellranger_count_b4aggr.output.aggr_input_csv
    elif "cellranger_vdj" in rules_available:
        return rules.cellranger_vdj_b4aggr.output.aggr_input_csv
    else:
        return rules.cellranger_multi_b4aggr.output.aggr_input_csv


# Define the function to determine output for cellranger_aggr
def output_dir_for_cellranger_aggr(wc):
    rules_available = list(rules.__dict__.keys())
    if config_aggr.get("aggr_outdir", None):
        return config_aggr["aggr_outdir"]
    elif "cellranger_count" in rules_available:
        return config_count["output_count"] + "/aggr_results"
    elif "cellranger_vdj" in rules_available:
        return config_vdj["output_vdj"] + "/aggr_results"
    else:
        return config_multi["output_multi"] + "/aggr_results"


# write a multi csv input file
def multi_csv():
    rules_available = list(rules.__dict__.keys())


# Define the function to determine input for cellranger_aggr
def input_for_cellranger_multi(wc):
    if config_multi.get("csv", None):
        return config_aggr["csv"]
    else:
        return multi_csv(rules)
