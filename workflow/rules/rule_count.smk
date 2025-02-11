# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers
import os, sys
import pandas as pd

count_outdir = os.path.join(fastq_outdirectory, "count_results")
donor_list = list(set(pd.read_csv(config_count["add_info_aggr"])["donor"].tolist()))


rule fastq_folders:
    input:
        os.path.join(fastq_outdirectory, "sample_to_fastq.json"),
    output:
        expand(os.path.join(fastq_outdirectory, "{sample}_fastq.csv"), sample=samples),
    shell:
        """
        python scripts/get_proxy_fastq_files.py {input}
        """


rule cellranger_count_b4aggr:
    output:
        aggr_input_csv=os.path.join(count_outdir, "aggregation_count.csv"),
    input:
        count_info_h5=expand(
            os.path.join(count_outdir, "{sample}_count", "outs", "molecule_info.h5"),
            sample=samples,
        ),
        FilteredBCmatricesHDF5=expand(
            os.path.join(
                count_outdir, "{sample}_count", "outs", "filtered_feature_bc_matrix.h5"
            ),
            sample=samples,
        ),
    log:
        os.path.join(count_outdir, "logs", "count_aggr.log"),
    params:
        countdir=count_outdir,
        additional_info_aggr=config_count["add_info_aggr"],
    container:
        "docker://litd/docker-cellranger:v8.0.1"
        # "cellranger.v8.0.1.sif"
    shell:
        """
        bash scripts/get_count_aggr_csv.sh {params.countdir} > {params.countdir}/aggregation_count.csv;
        python scripts/add_info_aggr.py {params.countdir}/aggregation_count.csv {params.additional_info_aggr}
        """


# Rule to count features for single library
rule cellranger_count:
    input:
        os.path.join(fastq_outdirectory, "{sample}_fastq.csv"),
    output:
        Permoleculereadinformation=os.path.join(
            count_outdir, "{sample}_count", "outs", "molecule_info.h5"
        ),
        FilteredBCmatricesHDF5=os.path.join(
            count_outdir, "{sample}_count", "outs", "filtered_feature_bc_matrix.h5"
        ),
    params:
        transcriptome=config_count["transcriptome"],
        id2use=lambda wc: f"pipestance_{wc.sample}",
        sample=lambda wc: f"{wc.sample}",
        count_outdir2use=lambda wc: os.path.join(count_outdir, f"{wc.sample}_count"),
        input_dir=lambda wc: os.path.join(fastq_outdirectory, f"{wc.sample}_fastq"),
        add_arguments=config_count["additional_arguments"],
    resources:
        cores=config_count["resources"]["localcores"],
        memory=config_count["resources"]["localmem"],
    container:
        "docker://litd/docker-cellranger:v8.0.1"
        # "cellranger.v8.0.1.sif"
    log:
        file=os.path.join(count_outdir, "logs", "count_{sample}.log"),
    benchmark:
        os.path.join(count_outdir, "benchmarks", "benchmarks_{sample}_count.csv")
    shell:
        """   
        cellranger count \
        --id={params.id2use} \
        --transcriptome={params.transcriptome} \
        --fastqs={params.input_dir} \
        --sample={params.sample} \
        --localcores={resources.cores} \
        --localmem={resources.memory} \
        {params.add_arguments} \
        2>&1 | tee -a {log.file}; \
        bash scripts/move_pipestance_count_dir.sh {log.file} {params.count_outdir2use}; 
        """
