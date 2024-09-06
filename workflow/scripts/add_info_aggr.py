import pandas as pd
import sys, os, re


def add_columns(org_csv, add_csv):

    org_data = pd.read_csv(org_csv)
    add_data = pd.read_csv(add_csv)

    new_data = pd.merge(org_data, add_data, how="left", on="sample_id")

    new_data.to_csv(org_csv, quotechar='"', index=False)


if __name__ == "__main__":

    original_aggr_csv = sys.argv[1]

    if len(sys.argv) == 3:
        additional_info_csv = sys.argv[2]
        add_columns(original_aggr_csv, additional_info_csv)
