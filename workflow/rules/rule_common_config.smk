import yaml
import os, sys, argparse
from scripts import create_multi_csv
import pandas as pd

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


def extract_data_block(file_to_extract):

    # Open the file and read lines
    df = pd.read_csv(file_to_extract)

    # Find the index of the "[Data]" block
    data_index = None
    for i, line in enumerate(df.iloc[:, 0].tolist()):
        if str(line).startswith("[Data]"):
            data_index = i
            break

    if data_index is None:
        raise ValueError("No '[Data]' block found in the file.")

    # Read rows after the "[Data]" block into a DataFrame
    data_lines = data_index + 2
    data = pd.read_csv(file_to_extract, skiprows=data_lines)

    return data


def create_sample_dict_from_csv(file_path):
    """
    Reads a CSV file and creates a dictionary with 'donor' as keys and
    lists of samples (e.g., sample_1, sample_2, ...) as values.

    Parameters:
        file_path (str): Path to the CSV file.

    Returns:
        dict: A dictionary with donors as keys and sample lists as values.
    """
    # Read the CSV file
    df = pd.read_csv(file_path, index_col=0)

    # Convert the DataFrame to a dictionary
    sample_dict = df.apply(lambda row: row.dropna().tolist(), axis=1).to_dict()

    return sample_dict
