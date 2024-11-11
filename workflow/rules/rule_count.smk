# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd

count_outdir = config_count["output_count"]

# fastq_df = pd.read_csv(config_count["paths2fastq_file"])
if config_count.get("paths2fastq_file") is not None:
    fastq_file_path = config_count["paths2fastq_file"]
else:
    bcl_run_index1 = [
        i
        for i, v in feature_type_dict.items()
        if v
        in [
            "Gene Expression",
            "gene expression",
            "gene-expression",
            "Gene-Expression",
        ]
    ]
    # print(bcl_run_index1)
    fastq_file_path = (
        fastq_outdirectory_dict[bcl_run_index1[0]]
        + "/mkfastq_success_"
        + bcl_run_index1[0]
        + ".csv"
    )
    # print(fastq_file_path)


fastq_df = pd.read_csv(fastq_file_path)
samples_seqs = fastq_df["sample"].tolist()
samples_paths = fastq_df["fastq"].tolist()
fastq_dict = {s: samples_paths[i] for i, s in enumerate(samples_seqs)}


rule cellranger_count_b4aggr:
    output:
        aggr_input_csv=expand(
            "{count_outdir}/aggregation_count.csv", count_outdir=count_outdir
        ),
    input:
        # flag=os.path.join(outdir, "mkfastq.sucess.csv"),
        count_info_h5=expand(
            "{count_outdir}/{sample}_count/outs/molecule_info.h5",
            sample=list(fastq_dict.keys()),
            count_outdir=count_outdir,
        ),
    params:
        countdir=count_outdir,
        additional_info_aggr=config_count["add_info_aggr"],
    shell:
        """
        bash scripts/get_count_aggr_csv.sh {params.countdir} > {params.countdir}/aggregation_count.csv;
        python scripts/add_info_aggr.py {params.countdir}/aggregation_count.csv {params.additional_info_aggr}
        """


# Rule to count features for single library
rule cellranger_count:
    input:
        fastq_folder=lambda wc: fastq_dict[wc.sample],
        transcriptome=config_count["transcriptome"],
    output:
        Permoleculereadinformation="{count_outdir}/{sample}_count/outs/molecule_info.h5",
        FilteredBCmatricesHDF5="{count_outdir}/{sample}_count/outs/filtered_feature_bc_matrix.h5",
    params:
        id2use="pipestance_{sample}",
        sample="{sample}",
        count_outdir2use="{count_outdir}/{sample}_count",
        add_arguments=config_count["additional_arguments"],
    resources:
        cores=config_count["resources"]["localcores"],
        memory=config_count["resources"]["localmem"],
    container:
        "docker://litd/docker-cellranger:v8.0.1"
    log:
        file="{count_outdir}/logs/count_{sample}.log",
    benchmark:
        "{count_outdir}/benchmarks/benchmarks_{sample}_count.csv"
    shell:
        """   
        cellranger count \
        --id={params.id2use} \
        --transcriptome={input.transcriptome} \
        --fastqs={input.fastq_folder} \
        --sample={params.sample} \
        --localcores={resources.cores} \
        --localmem={resources.memory} \
        {params.add_arguments} \
        2>&1 | tee -a {log.file}; \
        bash scripts/move_pipestance_count_dir.sh {log.file} {params.count_outdir2use}; 
        """
