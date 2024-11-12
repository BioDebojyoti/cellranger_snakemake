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
