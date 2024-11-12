import sys
import os
import re
import shutil
import pandas as pd


def generate_multi_input(config_multi, multi_outdir):

    if isinstance(config_multi.get("multi_config_csv"), list):

        config_files_for_multi = config_multi["multi_config_csv"]
        donor_list = []

        for curr_multi_file in config_files_for_multi:

            curr_donor = os.path.basename(curr_multi_file)
            curr_donor = re.sub("_multi_samplesheet.csv", "", curr_donor)
            donor_list.append(curr_donor)

            target_file = f"{curr_donor}_multi_samplesheet.csv"
            destination_path = os.path.join(multi_outdir, target_file)

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


if __name__ == "__main__":
    config_file = sys.argv[1]
    multi_outdir = sys.argv[2]

    generate_multi_input(config_multi=config_file, multi_outdir=multi_outdir)
