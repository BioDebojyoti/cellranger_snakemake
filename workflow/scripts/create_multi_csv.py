import csv
import argparse


def create_multi_csv(
    input_csvs=None, gene_expression_args=None, vdj_args=None, output_file="output.csv"
):
    # Prepare the basic structure for the CSV
    csv_content = []

    # Combine all rows from multiple CSVs
    libraries_data = []

    # print(input_csvs)

    # Process each input CSV file provided
    # if input_csvs:
    #     for input_csv in input_csvs:
    #         print(input_csv)
    #         with open(input_csv, "r") as f:
    #             reader = csv.DictReader(f)
    #             libraries_data.extend(list(reader))  # Read all rows and accumulate
    with open(input_csvs, "r") as f:
        reader = csv.DictReader(f)
        libraries_data.extend(list(reader))  # Read all rows and accumulate

    # Add Gene Expression block if gene_expression_args is provided
    if gene_expression_args:
        csv_content.append("[gene-expression],,")
        if gene_expression_args.get("reference"):
            csv_content.append(f"reference,\"{gene_expression_args['reference']}\",")
            csv_content.append(f"create-bam,{gene_expression_args['create_bam']},")
            csv_content.append(f",,")
        if gene_expression_args.get("probe_set"):
            csv_content.append(f"probe-set,{gene_expression_args['probe_set']},")
            csv_content.append(f",,")
        if gene_expression_args.get("filter_probes"):
            csv_content.append(
                f"filter-probes,{gene_expression_args['filter_probes']},"
            )
            csv_content.append(f",,")
        if gene_expression_args.get("r1_length"):
            csv_content.append(f"r1-length,{gene_expression_args['r1_length']},")
            csv_content.append(f",,")
        if gene_expression_args.get("r2_length"):
            csv_content.append(f"r2-length,{gene_expression_args['r2_length']},")
            csv_content.append(f",,")
        if gene_expression_args.get("chemistry"):
            csv_content.append(f"chemistry,{gene_expression_args['chemistry']},")
            csv_content.append(f",,")
        if gene_expression_args.get("expect_cells"):
            csv_content.append(f"expect-cells,{gene_expression_args['expect_cells']},")
            csv_content.append(f",,")
        if gene_expression_args.get("force_cells"):
            csv_content.append(f"force-cells,{gene_expression_args['force_cells']},")
            csv_content.append(f",,")
        if gene_expression_args.get("no_secondary"):
            csv_content.append(f"no-secondary,{gene_expression_args['no_secondary']},")
            csv_content.append(f",,")
        if gene_expression_args.get("check_compatibility"):
            csv_content.append(
                f"check-library-compatibility,{gene_expression_args['check_compatibility']},"
            )
            csv_content.append(f",,")
        if gene_expression_args.get("include_introns"):
            csv_content.append(
                f"include-introns,{gene_expression_args['include_introns']},"
            )
            csv_content.append(f",,")
        csv_content.append("")

    # Add VDJ block if vdj_args is provided
    if vdj_args:
        csv_content.append("[vdj]")
        if vdj_args.get("reference"):
            csv_content.append(f"reference,\"{vdj_args['reference']}\",")
            csv_content.append(f",,")
        if vdj_args.get("primers"):
            csv_content.append(f"inner-enrichment-primers,{vdj_args['primers']},")
            csv_content.append(f",,")
        if vdj_args.get("r1_length"):
            csv_content.append(f"r1-length,{vdj_args['r1_length']},")
            csv_content.append(f",,")
        if vdj_args.get("r2_length"):
            csv_content.append(f"r2-length,{vdj_args['r2_length']},")
            csv_content.append(f",,")
        csv_content.append("")

    # Add the Libraries block only if there are any libraries present
    if libraries_data:
        csv_content.append("[libraries],,")
        # Add the header line for the libraries block only once
        csv_content.append("fastq_id,fastqs,feature_types")

        for row in libraries_data:
            fastq_id = row["sample"]
            fastqs = row["fastq"]
            feature_type = row["lib_type"]
            csv_content.append(f"{fastq_id},{fastqs},{feature_type}")
        csv_content.append("")

    # print(output_file)
    # print(csv_content)
    # Write the final structure to the output file
    with open(output_file, "w") as f:
        f.write("\n".join(csv_content))


if __name__ == "__main__":
    # Set up the argument parser
    parser = argparse.ArgumentParser(
        description="Generate a multi-block CSV for libraries."
    )

    # Add argument for multiple input CSV files
    parser.add_argument(
        "--input-csvs",
        type=str,
        nargs="+",
        required=True,
        help="List of input CSV files.",
    )
    parser.add_argument(
        "--output-file",
        type=str,
        default="output.csv",
        help="Path to the output CSV file.",
    )

    # Gene Expression-specific arguments (optional)
    parser.add_argument(
        "--gex-reference",
        type=str,
        required=False,
        help="Path to the Gene Expression transcriptome reference file.",
    )
    parser.add_argument(
        "--gex-probe-set",
        type=str,
        required=False,
        help="Path to the probe set for Gene Expression.",
    )
    parser.add_argument(
        "--gex-filter-probes",
        type=str,
        required=False,
        help="Filter probes for Gene Expression (true/false).",
    )
    parser.add_argument(
        "--gex-r1-length",
        type=int,
        required=False,
        help="R1 length for Gene Expression.",
    )
    parser.add_argument(
        "--gex-r2-length",
        type=int,
        required=False,
        help="R2 length for Gene Expression.",
    )
    parser.add_argument(
        "--gex-chemistry",
        type=str,
        required=False,
        default="auto",
        help="Chemistry type for Gene Expression.",
    )
    parser.add_argument(
        "--gex-expect-cells",
        type=int,
        required=False,
        help="Expected number of cells for Gene Expression.",
    )
    parser.add_argument(
        "--gex-force-cells",
        type=int,
        required=False,
        help="Force cell count for Gene Expression.",
    )
    parser.add_argument(
        "--gex-no-secondary",
        type=str,
        required=False,
        help="No secondary alignment for Gene Expression (true/false).",
    )
    parser.add_argument(
        "--gex-create-bam",
        type=str,
        default=False,
        required=False,
        help="No BAM output for Gene Expression (true/false).",
    )
    parser.add_argument(
        "--gex-check-compatibility",
        type=str,
        required=False,
        help="Check library compatibility (true/false).",
    )
    parser.add_argument(
        "--gex-include-introns",
        type=str,
        required=False,
        help="Include introns (true/false).",
    )

    # VDJ-specific arguments (optional)
    parser.add_argument(
        "--vdj-reference",
        type=str,
        required=False,
        help="Path to the VDJ reference file.",
    )
    parser.add_argument(
        "--vdj-primers", type=str, required=False, help="Path to the VDJ primers file."
    )
    parser.add_argument(
        "--vdj-r1-length", type=int, required=False, help="R1 length for VDJ."
    )
    parser.add_argument(
        "--vdj-r2-length", type=int, required=False, help="R2 length for VDJ."
    )

    # Parse the arguments
    args = parser.parse_args()

    # Prepare argument dictionaries only if the respective arguments are provided
    gene_expression_args = (
        {
            "reference": args.gex_reference,
            "probe_set": args.gex_probe_set,
            "filter_probes": args.gex_filter_probes,
            "r1_length": args.gex_r1_length,
            "r2_length": args.gex_r2_length,
            "chemistry": args.gex_chemistry,
            "expect_cells": args.gex_expect_cells,
            "force_cells": args.gex_force_cells,
            "no_secondary": args.gex_no_secondary,
            "create_bam": args.gex_create_bam,
            "check_compatibility": args.gex_check_compatibility,
            "include_introns": args.gex_include_introns,
        }
        if args.gex_reference
        else None
    )

    vdj_args = (
        {
            "reference": args.vdj_reference,
            "primers": args.vdj_primers,
            "r1_length": args.vdj_r1_length,
            "r2_length": args.vdj_r2_length,
        }
        if args.vdj_reference
        else None
    )

    # Call the function to create the multi-block CSV
    create_multi_csv(
        input_csvs=args.input_csvs,
        gene_expression_args=gene_expression_args,
        vdj_args=vdj_args,
        output_file=args.output_file,
    )
