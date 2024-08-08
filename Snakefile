# Snakemake pipeline for scRNA-seq data using Cell Ranger and Seurat with Docker containers

import pandas as pd 

# Configuration
configfile: "config.yaml"

# Load samples from CSV
samples_df = pd.read_csv(config['sample_csv'])
samples = samples_df['sample'].tolist()

# Define lambda functions to generate paths based on config
bcl_folder = lambda wildcards: config["bcl_folder"] if wildcards else None
fastq_folder = lambda wildcards: config["flowcell_id"] + "/outs/fastq_path" if wildcards else config["fastq_path"]

# # Define rules
rule all:
    input:
        expand("{sample}/outs/filtered_feature_bc_matrix.h5", sample=config["samples"]),
        config["aggregation_id"] + "/outs/filtered_feature_bc_matrix.h5" if config.get("aggregate", False) else [],
        "reanalyzed/outs/filtered_feature_bc_matrix.h5" if config.get("reanalyze", False) else []


# Rule to generate FASTQ files if starting from BCL files
rule mkfastq:
    input:
        f"{bcl_folder(config["start_with_bcl"])}"
    output:
        f"{fastq_folder(config["start_with_bcl"])}"

    resources:
        cores = config["resources"]["localcores"],
        mem_mb = config["resources"]["localmem"] * 1024  # Convert GB to MB     

    container:
        "docker://biodebojyoti/crossplatformcellranger:8.0.1"
    shell:
        """
        cellranger mkfastq --run {input} --id={config[flowcell_id]} --output-dir={output} \
        --localcores={resources.cores} \
        --localmem={resources.mem_mb} 
        """

# Rule to count features for single library
rule cellranger_count:
    input:
        fastq_folder = f"{fastq_folder(config["start_with_bcl"])}",
        transcriptome = config["transcriptome"]
    output:
        matrix_h5 = "{sample}/outs/filtered_feature_bc_matrix.h5"
    params:
        sample = "{sample}",
        expect_cells = config.get("expect_cells", None),
        pass_expect_cells = config.get("pass_expect_cells", False),
        force_cells = config.get("force_cells", None),
        pass_force_cells = config.get("pass_force_cells", False),
        feature_ref = config.get("feature_ref", None),
        use_feature_ref = config.get("use_feature_ref", False),
        libraries_csv = config.get("libraries_csv", None),
        use_libraries = config.get("use_libraries", False)        

    resources:
        cores = config["resources"]["localcores"],
        mem_mb = config["resources"]["localmem"] * 1024  # Convert GB to MB        
    
    container:
        "docker://biodebojyoti/crossplatformcellranger:8.0.1"

    shell:
        """
        # Check if Feature Barcode parameters are required
        if [ "{params.use_feature_ref}" = "true" ] && [ -n "{params.feature_ref}" ]; then
        feature_ref_param="--feature-ref={params.feature_ref}"
        else
        feature_ref_param=""
        fi

        if [ "{params.use_libraries}" = "true" ] && [ -n "{params.libraries_csv}" ]; then
        libraries_param="--libraries={params.libraries_csv}"
        else
        libraries_param=""
        fi

        # For Gene Expression with Feature Barcode
        if [ "{params.use_feature_ref}" = "true" ] && [ -n "{params.feature_ref}" ]; then
        cellranger count --id={params.sample} \
            --transcriptome={input.transcriptome} \
            $libraries_param \
            $feature_ref_param \
            --create-bam={config[create_bam]} \
            --localcores={resources.cores} \
            --localmem={resources.mem_mb} \
        else
        # For Gene Expression without Feature Barcode
        echo cellranger count --id={params.sample} \
            --transcriptome={input.transcriptome} \
            --fastqs={input.fastq_folder} \
            --sample={params.sample} \
            --create-bam={config[create_bam]} \
            --localcores={resources.cores} \
            --localmem={resources.mem_mb} \
            $expect_cells_param \
            $force_cells_param
        fi
        """


# Rule to aggregate libraries (optional)
rule cellranger_aggr:
    input:
        csv = "aggregation.csv",
    output:
        matrix_h5 = config["aggregation_id"] + "/outs/filtered_feature_bc_matrix.h5"

    resources:
        cores = config["resources"]["localcores"],
        mem_mb = config["resources"]["localmem"] * 1024  # Convert GB to MB        

    container:
        "docker://biodebojyoti/crossplatformcellranger:8.0.1"
    
    shell:
        """
        cellranger aggr --id=config["aggregation_id"] --csv={input.csv} --normalize=mapped \
        --localcores={resources.cores} \
        --localmem={resources.mem_mb}
        """