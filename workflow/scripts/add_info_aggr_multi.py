import pandas as pd
import sys, os, re
import csv


def add_columns(org_csv, add_csv):

    org_data = pd.read_csv(org_csv)
    add_data = pd.read_csv(add_csv)

    add_data.drop(["sample_id"], inplace=True, axis=1)

    add_data = add_data.groupby("donor").agg(lambda x: x.unique()[0])
    add_data.reset_index(inplace=True)

    new_data = pd.merge(
        org_data,
        add_data,
        how="inner",
        left_on="sample_id",
        right_on="donor",
        validate="one_to_one",
    )

    with open(org_csv, mode="w", newline="") as file:
        # Write the header without quotes
        file.write(",".join(new_data.columns) + "\n")
    new_data.to_csv(org_csv, quoting=csv.QUOTE_ALL, mode="a", index=False, header=False)


if __name__ == "__main__":

    original_aggr_csv = sys.argv[1]

    if len(sys.argv) == 3:
        additional_info_csv = sys.argv[2]
        add_columns(original_aggr_csv, additional_info_csv)
