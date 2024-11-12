import pandas as pd
import argparse


def create_donor_dict(additional_info_csv, fastq_csv_files, multi_csv_directory):
    # Read additional_info_csv file
    additional_info_df = pd.read_csv(additional_info_csv)

    # Read and concatenate the fastq CSV files
    dfs = [pd.read_csv(f) for f in fastq_csv_files]
    fastq_df = pd.concat(dfs, ignore_index=True)

    # Merge based on 'sample_id' from additional_info_df and 'sample' from fastq_df
    merged_df = pd.merge(
        additional_info_df,
        fastq_df,
        left_on="sample_id",
        right_on="sample",
        how="inner",
        validate="one_to_many",
    )

    # Create a dictionary for each donor
    donor_dict = {}
    for donor, group_df in merged_df.groupby("donor"):
        donor_dict[donor] = group_df.to_dict(orient="records")
        donor_df = pd.DataFrame.from_dict(donor_dict[donor]).loc[
            :, ["fastq", "sample", "library_type"]
        ]
        output_path = multi_csv_directory + "/" + donor + "_multi.csv"

        # print(output_path)

        donor_df.to_csv(output_path, index=False)

    return list(donor_dict.keys())


def main():
    # Argument parser
    parser = argparse.ArgumentParser(
        description="Create a dictionary of donors by merging CSV files."
    )
    parser.add_argument(
        "--additional_info_csv",
        type=str,
        required=True,
        help="Path to the additional_info_csv file (e.g., 'additional_info_csv.csv')",
    )
    parser.add_argument(
        "--fastq_csv_files",
        type=str,
        nargs="+",
        required=True,
        help="the fastq path CSV files",
    )

    parser.add_argument(
        "--multi_csv_directory",
        type=str,
        required=True,
        help="Paths to the created CSV files",
    )

    # Parse arguments
    args = parser.parse_args()

    # Create donor dictionary
    donor_dict = create_donor_dict(
        args.additional_info_csv, args.fastq_csv_files, args.multi_csv_directory
    )


if __name__ == "__main__":
    main()
