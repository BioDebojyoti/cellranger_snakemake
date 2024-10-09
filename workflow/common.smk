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
